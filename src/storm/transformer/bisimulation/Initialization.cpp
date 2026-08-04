#include "storm/transformer/bisimulation/Initialization.h"

#include <algorithm>
#include <map>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/logic/AtomicExpressionFormula.h"
#include "storm/logic/AtomicLabelFormula.h"
#include "storm/logic/Formula.h"
#include "storm/logic/FormulaInformation.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

template<typename ValueType>
PreservationInformation Initialization<ValueType>::getPreservationInformation() const {
    storm::bisimulation::PreservationInformation information;
    // Add all labels that appear in all formulas
    for (auto const& f : formulas) {
        f->gatherReferencedRewardModels(information.preservedRewardModels);
        for (auto const& l : f->getAtomicLabelFormulas()) {
            information.preservedStateLabels.insert(l->getLabel());
        }
        for (auto const& e : f->getAtomicExpressionFormulas()) {
            information.preservedStateLabels.insert(e->toString());
        }
        STORM_LOG_ASSERT(std::all_of(information.preservedStateLabels.begin(), information.preservedStateLabels.end(),
                                     [this](std::string const& label) { return model.getStateLabeling().containsLabel(label); }),
                         "Formula " << *f << " uses a label that is not known in the model.");
        STORM_LOG_ASSERT(std::all_of(information.preservedRewardModels.begin(), information.preservedRewardModels.end(),
                                     [this](std::string const& rew) { return model.getRewardModels().contains(rew); }),
                         "Formula " << *f << " uses a reward model that is not known in the model.");
    }
    // if requested or if there are no formulas given (and the user hasn't explicitly set something else), also add all state labels and rewards
    if (options.preserveAllStateLabels.value_or(formulas.empty())) {
        for (auto const& label : model.getStateLabeling().getLabels()) {
            information.preservedStateLabels.insert(label);
        }
    }
    if (options.preserveAllRewards.value_or(formulas.empty())) {
        for (auto const& [name, rm] : model.getRewardModels()) {
            information.preservedRewardModels.insert(name);
        }
    }
    // if requested and available, preserve all choice labels
    if (options.preserveAllChoiceLabels.value_or(false) && model.hasChoiceLabeling()) {
        for (auto const& label : model.getChoiceLabeling().getLabels()) {
            information.preservedChoiceLabels.insert(label);
        }
    }
    return information;
}

template<typename ValueType>
Initialization<ValueType>::Initialization(storm::models::sparse::Model<ValueType> const& model, Options const& options,
                                          std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas)
    : model(model), options(options), formulas(formulas) {
    // sanity checks
    for (auto const& f : this->formulas) {
        auto const fInfo = f->info();
        if (fInfo.containsBoundedUntilFormula() || fInfo.containsNextFormula() || fInfo.containsCumulativeRewardFormula() || fInfo.containsDiscountFormula()) {
            STORM_LOG_THROW(options.bisimulationType == Options::BisimulationType::Strong, storm::exceptions::IllegalFunctionCallException,
                            "The formula " << *f << " is not known to be preserved by weak bisimulation.");
            // TODO: e.g. time-bounded reachability for CTMC should work, right?
        }
    }
    auto const preservationInformation = getPreservationInformation();

    // otherwise go through all relevant labels / rewards / model components and add them to the corresponding preserved annotations.
    for (auto const& label : preservationInformation.preservedStateLabels) {
        preservedStateAnnotations.booleans.push_back(std::cref(model.getStateLabeling().getStates(label)));
    }
    if (model.hasChoiceLabeling()) {
        for (auto const& label : preservationInformation.preservedChoiceLabels) {
            preservedChoiceAnnotations.booleans.push_back(std::cref(model.getChoiceLabeling().getChoices(label)));
        }
    }
    for (auto const& rewName : preservationInformation.preservedRewardModels) {
        auto const& rewardModel = model.getRewardModel(rewName);
        STORM_LOG_THROW(!rewardModel.hasTransitionRewards(), storm::exceptions::NotSupportedException,
                        "Bisimulation initialization does not handle transition rewards of reward model '" << rewName << "'.");
        if (rewardModel.hasStateRewards()) {
            preservedStateAnnotations.values.emplace_back(rewardModel.getStateRewardVector());
        }
        if (rewardModel.hasStateActionRewards()) {
            preservedChoiceAnnotations.values.emplace_back(rewardModel.getStateActionRewardVector());
        }
    }
    using enum storm::models::ModelType;
    if (model.isOfType(MarkovAutomaton)) {
        auto const& ma = model.template as<storm::models::sparse::MarkovAutomaton<ValueType>>();
        preservedStateAnnotations.booleans.emplace_back(ma->getMarkovianStates());
        preservedStateAnnotations.values.emplace_back(ma->getExitRates());
    } else {
        STORM_LOG_THROW(model.isOfType(Dtmc) || model.isOfType(Ctmc) || model.isOfType(Mdp), storm::exceptions::NotSupportedException,
                        "Bisimulation initialization is not implemented for model type '" << model.getType() << "'.");
    }
}
template<typename ValueType>
bool Initialization<ValueType>::LocalSignatureComp::operator()(uint64_t a, uint64_t b) const {
    // handle action sensitive case
    if (!extraAnnotation.empty() && extraAnnotation[a] != extraAnnotation[b]) {
        return extraAnnotation[a] < extraAnnotation[b];
    }
    for (auto const& v : annotations.booleans) {
        if (bool v2 = v.get().get(b); v.get().get(a) != v2)
            return v2;
    }
    for (auto const& v : annotations.integers) {
        if (v[a] != v[b])
            return v[a] < v[b];
    }
    for (auto const& v : annotations.values) {
        if (v[a] != v[b])  // todo: add some floating point slack here?
            return v[a] < v[b];
    }
    return false;
}

template<typename ValueType>
std::vector<uint64_t> Initialization<ValueType>::getChoiceClasses() const {
    auto const& transitionMatrix = model.getTransitionMatrix();

    std::vector<uint64_t> choiceClasses;
    choiceClasses.reserve(model.getNumberOfChoices());

    // maintain a map from choice index to classes that distinguishes choices iff their preserved annotations are different.

    std::vector<uint64_t> localActionIndices;
    if (options.actionSensitive && !transitionMatrix.hasTrivialRowGrouping()) {
        localActionIndices.reserve(model.getNumberOfChoices());
        for (uint64_t state = 0; state < model.getNumberOfChoices(); ++state) {
            for (uint64_t act = 0; act < model.getTransitionMatrix().getRowGroupSize(state); ++act) {
                localActionIndices.push_back(act);
            }
        }
    }
    LocalSignatureComp comp(preservedChoiceAnnotations, localActionIndices);
    std::map<uint64_t, uint64_t, LocalSignatureComp> classes(comp);
    for (uint64_t state = 0; state < model.getNumberOfStates(); ++state) {
        for (auto const choice : transitionMatrix.getRowGroupIndices(state)) {
            // if we find a new class of choice, add it to the map and increment the class counter.
            auto const [it, inserted] = classes.try_emplace(choice, classes.size());
            choiceClasses.push_back(it->second);
        }
    }
    return choiceClasses;
}

template<typename ValueType>
std::vector<uint64_t> Initialization<ValueType>::getStateClasses(std::vector<uint64_t> const& choiceClasses) const {
    auto const& rowGroupIndices = model.getTransitionMatrix().getRowGroupIndices();
    uint64_t const numStates = model.getNumberOfStates();

    // We do this in two steps. First, we distinguish the states only by their (sorted) choice classes
    std::vector<uint64_t> choiceBasedClasses;
    {
        choiceBasedClasses.reserve(model.getNumberOfStates());
        std::map<std::vector<uint64_t>, uint64_t> cbClasses;
        for (uint64_t state = 0; state < numStates; ++state) {
            std::vector<uint64_t> orderedChoiceClasses(choiceClasses.begin() + rowGroupIndices[state], choiceClasses.begin() + rowGroupIndices[state + 1]);
            std::sort(orderedChoiceClasses.begin(), orderedChoiceClasses.end());
            auto const [it, inserted] = cbClasses.try_emplace(std::move(orderedChoiceClasses), cbClasses.size());
            choiceBasedClasses.push_back(it->second);
        }
    }

    std::vector<uint64_t> stateClasses;
    stateClasses.reserve(numStates);
    LocalSignatureComp comp(preservedStateAnnotations, choiceBasedClasses);
    std::map<uint64_t, uint64_t, LocalSignatureComp> classes(comp);
    for (uint64_t state = 0; state < model.getNumberOfStates(); ++state) {
        auto const [it, inserted] = classes.try_emplace(state, classes.size());
        stateClasses.push_back(it->second);
    }
    return stateClasses;
}

template class Initialization<double>;
template class Initialization<storm::RationalNumber>;
template class Initialization<storm::RationalFunction>;
template class Initialization<storm::Interval>;
template class Initialization<storm::RationalInterval>;

}  // namespace storm::bisimulation

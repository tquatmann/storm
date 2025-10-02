#pragma once

#include <memory>
#include <string>
#include <vector>

#include "storm/modelchecker/helper/ltl/internal/RabinObjective.h"
#include "storm/models/sparse/StateLabeling.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

#include "storm/exceptions/NotImplementedException.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"

namespace storm {

namespace transformer {
namespace detail {

using RabinCondition = storm::modelchecker::helper::RabinCondition;
using RabinObjective = storm::modelchecker::helper::RabinObjective;

struct AccMec {
    storm::storage::MaximalEndComponent component;
    storm::storage::BitVector obj;
};

template<typename ValueType>
struct AcceptingMecsCollector {
    storm::storage::SparseMatrix<ValueType> const &transitionMatrix, backwardTransitions;
    storm::models::sparse::StateLabeling const& stateLabeling;
    std::vector<RabinObjective> const &qualitativeObjectives, quantitativeObjectives;

    // for each found accepting MEC, a representative state is chosen. maps each state to the accepting MECs it represents.
    std::map<uint64_t, std::vector<AccMec>> stateToMACs;

    /*!
     * Check if the inf condition of the provided Rabin condition is satisfied in the given MEC, i.e., whether all sets are visited infinitely often.
     */
    bool satisfiesInfCondition(storm::storage::MaximalEndComponent const& ec, RabinCondition const& condition) const {
        return std::all_of(condition.inf.begin(), condition.inf.end(),
                           [&ec, this](auto const& label) { return ec.containsAnyState(stateLabeling.getStates(label)); });
    }

    void gatherAcceptingSubMecs(storm::storage::MaximalEndComponent const& mec, storm::storage::BitVector enabledqualitativeObjectives = {},
                                uint64_t objIndex = 0, std::function<bool(storm::storage::MaximalEndComponent const&)> const& infCallback = {}) {
        STORM_LOG_ASSERT(!infCallback || infCallback(mec), "MEC does not satisfy the infCallback");
        enabledqualitativeObjectives.resize(quantitativeObjectives.size(), false);

        if (objIndex >= qualitativeObjectives.size() + quantitativeObjectives.size()) {
            // all objectives processed, store the MEC
            STORM_LOG_ASSERT(enabledqualitativeObjectives.size() == quantitativeObjectives.size(), "unexpected size of enabled qualitative objectives");
            auto const representativeState = mec.begin()->first;
            stateToMACs[representativeState].push_back(AccMec{mec, std::move(enabledqualitativeObjectives)});
            return;
        }

        bool const isQualitative = objIndex < qualitativeObjectives.size();
        auto const& objective = isQualitative ? qualitativeObjectives[objIndex] : quantitativeObjectives[objIndex - qualitativeObjectives.size()];
        if (!isQualitative) {
            enabledqualitativeObjectives.set(objIndex - qualitativeObjectives.size(), true);
        }

        // Handle buchi and non-buchi conditions separately
        // The former can be handled without modifying the mec, i.e., we only need a single recursive call for them.
        std::vector<std::reference_wrapper<RabinCondition const>> satisfiableBuchiConditions;  // conditions without a fin-state in this mec
        for (auto const& cond : objective) {
            if (!satisfiesInfCondition(mec, cond)) {
                continue;  // the condition cannot be satisfied in this mec and will be ignored from now on.
            }
            if (cond.fin.has_value() && mec.containsAnyState(stateLabeling.getStates(cond.fin.value()))) {
                // we have a non-buchi rabin condition.
                // Drop the states that we shall only visit finitely often and recurse on the resulting accepting sub-mecs
                auto const& finStates = stateLabeling.getStates(cond.fin.value());
                storm::storage::BitVector allowedMecStates(transitionMatrix.getRowGroupCount(), false);
                for (auto const& [state, _] : mec) {
                    if (!finStates.get(state)) {
                        allowedMecStates.set(state, true);
                    }
                }
                auto subMecs = storm::storage::MaximalEndComponentDecomposition(transitionMatrix, backwardTransitions, allowedMecStates);
                for (auto const& subMec : subMecs) {
                    if (satisfiesInfCondition(subMec, cond)) {
                        // The sub-mec is accepting for this objective, so we will recurse on it
                        auto const subMecCallback = [&infCallback, &cond, this](auto const& m) {
                            return (!infCallback || infCallback(m)) && satisfiesInfCondition(m, cond);
                        };
                        gatherAcceptingSubMecs(subMec, enabledqualitativeObjectives, objIndex + 1, subMecCallback);
                    }
                }
            } else {
                // we have a Buchi condition which we handle below with a single recursive call
                satisfiableBuchiConditions.push_back(cond);
            }
        }

        if (!satisfiableBuchiConditions.empty()) {
            // gather accepting MECs with one of the buchi conditions.
            auto acceptedWithBuchi = [&infCallback, &satisfiableBuchiConditions, this](auto const& subMec) {
                return (!infCallback || infCallback(subMec)) &&
                       std::any_of(satisfiableBuchiConditions.begin(), satisfiableBuchiConditions.end(),
                                   [&subMec, this](auto const& buchi) { return satisfiesInfCondition(subMec, buchi.get()); });
            };
            gatherAcceptingSubMecs(mec, enabledqualitativeObjectives, objIndex + 1, acceptedWithBuchi);
        }

        if (!isQualitative) {
            // recursive call without caring for this objective
            enabledqualitativeObjectives.set(objIndex - qualitativeObjectives.size(), false);
            gatherAcceptingSubMecs(mec, enabledqualitativeObjectives, objIndex + 1, infCallback);
        }
    }
};

/*!
 * For each MAC, a copy of the MAC is created and a choice from the original representative state to the corresponding copy state is added.
 * The MACs are appended in the sense that all state indices and local ("offset") choice indices of the input model remain valid.
 * @tparam SparseModelType
 * @param inputModel
 * @param stateToMACs
 * @return
 */
template<typename SparseModelType>
std::shared_ptr<SparseModelType> appendMacCopies(SparseModelType const& inputModel, std::map<uint64_t, std::vector<AccMec>> const& stateToMACs,
                                                 std::vector<std::string> const& quantitativeTotalRewardModelNames) {
    STORM_LOG_THROW(inputModel.isOfType(storm::models::ModelType::Mdp), storm::exceptions::NotImplementedException, "Functionality only implemented for MDPs.");
    // Note: Markov Automata are not supported as the added choices to enter the MAC copies are not allowed on Markovian representative states.
    using ValueType = typename SparseModelType::ValueType;
    storm::storage::SparseMatrix<ValueType> const& inputTransitions = inputModel.getTransitionMatrix();

    // We first construct
    // (a) the transition matrix of the new model,
    // (b) action reward vectors that assign a reward of 1 to the choices that enter a MAC, and
    // (c) mappings to allow lifting other model components
    std::vector<std::vector<ValueType>> rabinRewards(quantitativeTotalRewardModelNames.size());  // resized below
    std::vector<uint64_t> toInputStateMap, toInputChoiceMap;
    uint64_t const InvalidInputChoice = std::numeric_limits<uint64_t>::max();
    storm::storage::SparseMatrixBuilder<ValueType> builder(0, 0, 0, true, true, 0);

    // add the transitions of the input model
    uint64_t firstStateOfNextMacCopy = inputTransitions.getRowGroupCount();
    uint64_t newChoice{0};
    uint64_t state{0};
    for (; state < inputTransitions.getRowGroupCount(); ++state) {
        toInputStateMap.push_back(state);
        builder.newRowGroup(newChoice);
        for (auto const& inputChoice : inputTransitions.getRowGroupIndices(state)) {
            toInputChoiceMap.push_back(inputChoice);
            for (auto const& entry : inputTransitions.getRow(inputChoice)) {
                builder.addNextValue(newChoice, entry.getColumn(), entry.getValue());
            }
            ++newChoice;
        }
        if (auto findRes = stateToMACs.find(state); findRes != stateToMACs.end()) {
            // add choice to the copy of the MAC
            for (auto const& mac : findRes->second) {
                toInputChoiceMap.push_back(InvalidInputChoice);
                builder.addNextValue(newChoice, firstStateOfNextMacCopy, storm::utility::one<ValueType>());
                for (auto const accRabinObjIndex : mac.obj) {  // assign reward 1 to all objectives that are accepting in this MAC
                    rabinRewards[accRabinObjIndex].resize(newChoice, storm::utility::zero<ValueType>());
                    rabinRewards[accRabinObjIndex].push_back(storm::utility::one<ValueType>());
                }
                ++newChoice;
                firstStateOfNextMacCopy += mac.component.size();
            }
        }
    }
    // add the MAC copies
    // This code assumes that the macs are enumerated in the same order as the choices to enter them above.
    // This is the case as stateToMACs is a map and thus the states are enumerated in ascending order.
    // allocate storage for a mapping of column indices of the input model to the corresponding column indices of the new model
    std::vector<uint64_t> inputToNewStateMap(inputTransitions.getColumnCount(), std::numeric_limits<uint64_t>::max());
    for (auto const& [reprState, macs] : stateToMACs) {
        for (auto const& mac : macs) {
            for (auto const& [macState, macChoices] : mac.component) {
                inputToNewStateMap[macState] = state;
                ++state;
            }
            for (auto const& [macState, macChoices] : mac.component) {
                toInputStateMap.push_back(macState);
                builder.newRowGroup(newChoice);
                for (auto const& macChoice : macChoices) {
                    toInputChoiceMap.push_back(macChoice);
                    for (auto const& entry : inputTransitions.getRow(macChoice)) {
                        builder.addNextValue(newChoice, inputToNewStateMap[entry.getColumn()], entry.getValue());
                    }
                    ++newChoice;
                }
            }
        }
    }
    STORM_LOG_ASSERT(state == firstStateOfNextMacCopy, "Inconsistent state count.");
    STORM_LOG_ASSERT(state == toInputStateMap.size(), "Inconsistent state count.");
    STORM_LOG_ASSERT(newChoice == toInputChoiceMap.size(), "Inconsistent choice count.");

    // helper lambda to lift vectors from the input model to the new model
    auto liftVector = [](auto const& inputVec, std::vector<uint64_t> const& toInputMap) {
        using VecType = std::decay_t<decltype(inputVec)>;
        auto const size = toInputMap.size();
        VecType result;
        if constexpr (std::is_same_v<VecType, storm::storage::BitVector>) {
            result.resize(size, false);
            for (uint64_t i = 0; i < size; ++i) {
                if (toInputMap[i] < inputVec.size() && inputVec.get(toInputMap[i])) {
                    result.set(i, true);
                }
            }
        } else {
            result.reserve(size);
            for (uint64_t i = 0; i < size; ++i) {
                if (toInputMap[i] < inputVec.size()) {
                    result.push_back(inputVec[toInputMap[i]]);
                } else {
                    result.push_back(storm::utility::zero<typename VecType::value_type>());
                }
            }
        }
        return result;
    };

    // Initialize the new model components
    uint64_t const stateCount = toInputStateMap.size();
    uint64_t const choiceCount = toInputChoiceMap.size();
    storm::storage::sparse::ModelComponents<ValueType> components(builder.build(choiceCount, stateCount, stateCount),
                                                                  storm::models::sparse::StateLabeling(state));
    // State labeling
    for (auto const& label : inputModel.getStateLabeling().getLabels()) {
        auto const& states = inputModel.getStateLabeling().getStates(label);
        if (label == "init") {
            auto initialStates = states;
            initialStates.resize(stateCount, false);  // init label not applied to copies of the MACs
            components.stateLabeling.addLabel(label, std::move(initialStates));
        } else {
            components.stateLabeling.addLabel(label, liftVector(states, toInputStateMap));
        }
    }
    // Reward models
    for (auto const& [rewName, rewModel] : inputModel.getRewardModels()) {
        std::optional<std::vector<ValueType>> newStateRew, newStateActionRew;
        if (rewModel.hasStateRewards()) {
            newStateRew = liftVector(rewModel.getStateRewardVector(), toInputStateMap);
        }
        if (rewModel.hasStateActionRewards()) {
            newStateActionRew = liftVector(rewModel.getStateActionRewardVector(), toInputChoiceMap);
        }
        STORM_LOG_THROW(!rewModel.hasTransitionRewards(), storm::exceptions::NotSupportedException, "Transition rewards are not supported.");
        storm::models::sparse::StandardRewardModel<ValueType> productRewModel(std::move(newStateRew), std::move(newStateActionRew));
        components.rewardModels.emplace(rewName, std::move(productRewModel));
    }
    for (uint64_t i = 0; i < quantitativeTotalRewardModelNames.size(); ++i) {
        rabinRewards[i].resize(choiceCount, storm::utility::zero<ValueType>());
        storm::models::sparse::StandardRewardModel<ValueType> rabinRewardModel(std::nullopt, std::move(rabinRewards[i]));
        auto [_, inserted] = components.rewardModels.emplace(quantitativeTotalRewardModelNames[i], std::move(rabinRewardModel));
        STORM_LOG_THROW(inserted, storm::exceptions::NotImplementedException, "Reward model name clash: " << quantitativeTotalRewardModelNames[i]);
    }
    return storm::utility::builder::buildModelFromComponents(inputModel.getType(), std::move(components));
}

}  // namespace detail

/*!
 * Demerges MECs
 */
class RabinToTotalRewardTransformer {
   public:
    using RabinCondition = detail::RabinCondition;
    using RabinObjective = detail::RabinObjective;

    template<typename SparseModelType>
    struct ReturnType {
        std::shared_ptr<SparseModelType> model;
        std::string finStatesLabel;
    };

    /*!
     * Produces a model with total reward objectives as follows.
     * First, a set of maximal accepting end components (MACs) is found such that
     * - All MACs satisfy the rabinObjectives for which index the qualitativeRabinObjectives parameter is true (assuming that each state in that MAC is
     * visited infinitely often).
     * - and a subset of the quantitative rabin objectives can be satisfied.
     * Then, for each MAC a copy of the MAC is created and a choice from one representative state of the MAC to the corresponding copy state is added.
     * Entering a MAC for such a choice gives a reward of 1 for each quantitative RabinObjective that is satisfied when each MAC state is visited infinitely
     * often. Here, the i'th total reward model name corresponds to the i'th quantitativeRabinObjective, i.e., the i'th 'false' entry in the
     * qualitativeRabinObjectives parameter.
     * The returned finStatesLabel indicate those states where it is sufficient (and in the case of qualitative objectives
     * also necessary) to visit them only finitely often.
     *
     * @return
     */
    template<typename SparseModelType>
    static ReturnType<SparseModelType> transform(SparseModelType const& model, std::vector<RabinObjective> const& rabinObjectives,
                                                 std::vector<std::string> const& quantitativeTotalRewardModelNames,
                                                 storm::storage::BitVector const& qualitativeRabinObjectives) {
        STORM_LOG_ASSERT(model.isNondeterministicModel(), "Expecting a nondeterministic model.");
        STORM_LOG_ASSERT(rabinObjectives.size() == qualitativeRabinObjectives.size(), "Inconsistent objective counts.");
        STORM_LOG_ASSERT(rabinObjectives.size() - qualitativeRabinObjectives.getNumberOfSetBits() == quantitativeTotalRewardModelNames.size(),
                         "Inconsistent number of quantitative rabin objectives.");
        auto qualitative = storm::utility::vector::filterVector(rabinObjectives, qualitativeRabinObjectives);
        auto quantitative = storm::utility::vector::filterVector(rabinObjectives, ~qualitativeRabinObjectives);
        auto backwardTransitions = model.getBackwardTransitions();

        // Identify MECs
        storm::storage::MaximalEndComponentDecomposition mecs(model.getTransitionMatrix(), backwardTransitions);

        // gather accepting (sub-)MECs
        detail::AcceptingMecsCollector accMecs{model.getTransitionMatrix(), backwardTransitions, model.getStateLabeling(), qualitative, quantitative};
        for (auto& mec : mecs) {
            accMecs.gatherAcceptingSubMecs(mec);
        }

        // build the new model
        ReturnType<SparseModelType> result;
        result.model = detail::appendMacCopies(model, accMecs.stateToMACs, quantitativeTotalRewardModelNames);

        // label the states that we shall visit only finitely often
        storm::storage::BitVector finStates(model.getNumberOfStates(), true);
        finStates.resize(result.model->getNumberOfStates(), false);
        result.finStatesLabel = result.model->getStateLabeling().addUniqueLabel("fin", std::move(finStates));
        return result;
    }
};
}  // namespace transformer
}  // namespace storm
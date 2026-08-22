#include "storm/transformer/bisimulation/Quotient.h"

#include <algorithm>
#include <limits>
#include <map>
#include <optional>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/sparse/ModelComponents.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/utility/builder.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

template<QuotientType quotientType, typename ValueType>
typename Quotient<quotientType, ValueType>::IndexMapping Quotient<quotientType, ValueType>::computeIndexMappings(
    storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition const& partition,
    std::optional<std::vector<uint64_t>> const& choiceClasses) {
    uint64_t const undef = std::numeric_limits<uint64_t>::max();
    IndexMapping indexMapping;
    indexMapping.toQuotientState.assign(model.getNumberOfStates(), undef);
    indexMapping.toRepresentativeState.reserve(partition.getNumberOfBlocks());
    if (model.isNondeterministicModel()) {
        indexMapping.toRepresentativeChoice.emplace();
    }
    uint64_t quotientState = 0;
    partition.forEachBlock([&](auto const& block) {
        for (auto const s : block) {
            indexMapping.toQuotientState[s] = quotientState;
        }
        indexMapping.toRepresentativeState.push_back(block.front());
        if (model.isNondeterministicModel()) {
            // TODO: Deduplicate based on choiceClasses and signature?
            for (auto const representativeChoice : model.getTransitionMatrix().getRowGroupIndices(block.front())) {
                indexMapping.toRepresentativeChoice->push_back(representativeChoice);
            }
        }
        ++quotientState;
    });
    STORM_LOG_ASSERT(std::none_of(indexMapping.toQuotientState.begin(), indexMapping.toQuotientState.end(), [&undef](auto const& s) { return s == undef; }),
                     "Not all states appear in a block of the partition.");
    if (model.isNondeterministicModel()) {
        indexMapping.toRepresentativeChoice->shrink_to_fit();
    }
    return indexMapping;
}

template<QuotientType quotientType, typename ValueType>
auto Quotient<quotientType, ValueType>::buildFromPartition(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Options const& options,
                                                           storm::bisimulation::PreservationInformation const& preservationInformation,
                                                           IndexMapping const& indexMapping)
    -> std::shared_ptr<storm::models::sparse::Model<QuotientValueType>> {
    STORM_LOG_THROW(options.bisimulationType == Options::BisimulationType::Strong, storm::exceptions::NotImplementedException,
                    "Weak bisimulation is not implemented.");
    STORM_LOG_ASSERT(!model.isNondeterministicModel() || indexMapping.toRepresentativeChoice.has_value(), "Empty choice mapping for nondeterministic model.");
    auto const& toQuotientState = indexMapping.toQuotientState;
    auto const& toRepresentativeState = indexMapping.toRepresentativeState;
    auto const& toRepresentativeChoice =
        indexMapping.toRepresentativeChoice.has_value() ? indexMapping.toRepresentativeChoice.value() : indexMapping.toRepresentativeState;

    uint64_t const numberOfQuotientStates = toRepresentativeState.size();
    uint64_t const numberOfQuotientChoices = toRepresentativeChoice.size();
    bool const isNondeterministic = model.isNondeterministicModel();
    STORM_LOG_ASSERT(isNondeterministic || numberOfQuotientStates == numberOfQuotientChoices, "Unexpected choice count.");

    // Now build the model components one after the other
    storm::storage::sparse::ModelComponents<QuotientValueType> components;

    // Build the transition matrix
    {
        storm::storage::SparseMatrixBuilder<QuotientValueType> builder(numberOfQuotientChoices, numberOfQuotientStates, 0, true, isNondeterministic,
                                                                       isNondeterministic ? numberOfQuotientStates : 0);
        for (uint64_t quotientState = 0, quotientChoice = 0; quotientState < numberOfQuotientStates; ++quotientState) {
            if (isNondeterministic) {
                builder.newRowGroup(quotientChoice);
            }
            // TODO: handle deduplication
            for (auto const representativeChoice : model.getTransitionMatrix().getRowGroupIndices(toRepresentativeState[quotientState])) {
                std::map<uint64_t, QuotientValueType> quotientRow;
                for (auto const& entry : model.getTransitionMatrix().getRow(representativeChoice)) {
                    if (auto const ret = quotientRow.emplace(toQuotientState[entry.getColumn()], entry.getValue()); !ret.second) {
                        ret.first->second += entry.getValue();
                    }
                }
                for (auto const& [column, value] : quotientRow) {
                    builder.addNextValue(quotientChoice, column, value);
                }
                ++quotientChoice;
            }
        }
        components.transitionMatrix = builder.build();
    }

    // build state labeling
    {
        components.stateLabeling = storm::models::sparse::StateLabeling(numberOfQuotientStates);
        for (auto const& l : preservationInformation.preservedStateLabels) {
            auto const& in = model.getStateLabeling().getStates(l);
            storm::storage::BitVector out(numberOfQuotientStates, false);
            for (uint64_t quotientState = 0; quotientState < numberOfQuotientStates; ++quotientState) {
                if (in.get(toRepresentativeState[quotientState])) {
                    out.set(quotientState, true);
                }
            }
            components.stateLabeling.addLabel(l, std::move(out));
        }
        if (!preservationInformation.preservedStateLabels.contains("init")) {
            storm::storage::BitVector init(numberOfQuotientStates, false);
            for (auto const i : model.getStateLabeling().getStates("init")) {
                init.set(toQuotientState[i], true);
            }
            components.stateLabeling.addLabel("init", std::move(init));
        }
    }

    // build choice labeling
    if (!preservationInformation.preservedChoiceLabels.empty()) {
        STORM_LOG_ASSERT(model.hasChoiceLabeling(), "Model has no choice labeling but bisimulation preserved some.");
        components.choiceLabeling.emplace(numberOfQuotientChoices);
        for (auto const& l : preservationInformation.preservedChoiceLabels) {
            auto const& in = model.getChoiceLabeling().getChoices(l);
            storm::storage::BitVector out(numberOfQuotientChoices, false);
            for (uint64_t quotientChoice = 0; quotientChoice < numberOfQuotientChoices; ++quotientChoice) {
                if (in.get(toRepresentativeChoice[quotientChoice])) {
                    out.set(quotientChoice, true);
                }
            }
            components.choiceLabeling->addLabel(l, std::move(out));
        }
    }

    // build choice origins
    if (model.hasChoiceOrigins() && options.preserveChoiceOrigins) {
        components.choiceOrigins = model.getChoiceOrigins()->selectChoices(toRepresentativeChoice);
    }

    // build reward models
    {
        for (auto const& r : preservationInformation.preservedRewardModels) {
            std::optional<std::vector<QuotientValueType>> stateRewards, stateActionRewards;
            auto const& rm = model.getRewardModel(r);
            if (rm.hasStateRewards()) {
                stateRewards.emplace();
                stateRewards->reserve(numberOfQuotientStates);
                for (auto const representativeState : toRepresentativeState) {
                    stateRewards->push_back(rm.getStateReward(representativeState));
                }
            }
            if (rm.hasStateActionRewards()) {
                stateActionRewards.emplace();
                stateActionRewards->reserve(numberOfQuotientChoices);
                for (auto const representativeChoice : toRepresentativeChoice) {
                    stateActionRewards->push_back(rm.getStateActionReward(representativeChoice));
                }
            }
            STORM_LOG_THROW(!rm.hasTransitionRewards(), storm::exceptions::NotSupportedException,
                            "Transition rewards are not supported for quotient construction");
            components.rewardModels.emplace(
                r, storm::models::sparse::StandardRewardModel<QuotientValueType>(std::move(stateRewards), std::move(stateActionRewards)));
        }
    }

    // build model type specific components
    using enum storm::models::ModelType;
    if (model.isOfType(Ctmc)) {
        components.rateTransitions = true;
    } else if (model.isOfType(MarkovAutomaton)) {
        auto const& ma = model.template as<storm::models::sparse::MarkovAutomaton<ValueType>>();
        components.markovianStates.emplace(numberOfQuotientStates, false);
        components.exitRates.emplace(numberOfQuotientStates, storm::utility::zero<QuotientValueType>());
        for (uint64_t quotientState = 0; quotientState < numberOfQuotientStates; ++quotientState) {
            auto const representativeState = toRepresentativeState[quotientState];
            if (ma->isMarkovianState(representativeState)) {
                components.markovianStates->set(quotientState, true);
                STORM_LOG_ASSERT(components.transitionMatrix.getRowGroupSize(quotientState) == 1, "Unexpected number of choices for Markovian state.");
                // TODO: At this point, we ignore potential floating point inaccuracies and just pick the rate of some representative state.
                // We could at least pick a middle value from the interval of all exit rates.
                // Update: The initial partition should have already taken care of only grouping states with approximately(?) the same exit rate.
                components.exitRates.value()[quotientState] = ma->getExitRate(representativeState);
            }
        }
    } else {
        STORM_LOG_THROW(model.isOfType(Dtmc) || model.isOfType(Mdp), storm::exceptions::NotSupportedException,
                        "The model type " << model.getType() << " is not supported for quotient construction.");
    }
    return storm::utility::builder::buildModelFromComponents(model.getType(), std::move(components));
}

template class Quotient<QuotientType::Exact, double>;
template class Quotient<QuotientType::Exact, storm::RationalNumber>;
template class Quotient<QuotientType::Exact, storm::RationalFunction>;
template class Quotient<QuotientType::Exact, storm::Interval>;
template class Quotient<QuotientType::Exact, storm::RationalInterval>;
template class Quotient<QuotientType::IntervalAbstraction, double>;
template class Quotient<QuotientType::IntervalAbstraction, storm::RationalNumber>;
template class Quotient<QuotientType::IntervalAbstraction, storm::Interval>;
template class Quotient<QuotientType::IntervalAbstraction, storm::RationalInterval>;

}  // namespace storm::bisimulation

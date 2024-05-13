#include "storm-pomdp/beliefs/exploration/BeliefExploration.h"

#include "storm-pomdp/beliefs/exploration/BeliefExplorationMatrix.h"
#include "storm-pomdp/beliefs/exploration/BeliefExplorationMode.h"
#include "storm-pomdp/beliefs/exploration/ExplorationQueue.h"

#include "storm-pomdp/beliefs/abstraction/FreudenthalTriangulationBeliefAbstraction.h"
#include "storm-pomdp/beliefs/abstraction/RewardBoundedBeliefSplitter.h"
#include "storm-pomdp/beliefs/exploration/FirstStateNextStateGenerator.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Pomdp.h"
#include "storm/utility/OptionalRef.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

namespace storm::pomdp::beliefs {

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
template<typename InfoType, typename NextStateHandleType>
bool BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::performExploration(InfoType& info, NextStateHandleType&& exploreNextStates,
                                                                                      TerminalBeliefCallback const& terminalBeliefCallback,
                                                                                      TerminationCallback const& terminationCallback) {
    while (info.queue.hasNext()) {
        // Check if we terminate prematurely
        if ((terminationCallback && terminationCallback()) || storm::utility::resources::isTerminate()) {
            return false;  // Terminate prematurely
        }

        // Get the next belief to explore and perform some checks
        auto const currentBeliefId = info.queue.popNext();
        STORM_LOG_ASSERT(info.discoveredBeliefs.containsId(currentBeliefId), "Unknown belief id");
        STORM_LOG_ASSERT(info.exploredBeliefs.count(currentBeliefId) == 0, "Belief #" << currentBeliefId << " already explored.");
        STORM_LOG_ASSERT(info.terminalBeliefValues.count(currentBeliefId) == 0, "Belief #" << currentBeliefId << " already found to be terminal.");
        // do not take the current belief as reference since it will be invalidated when collecting more beliefs
        auto const currentBelief = info.discoveredBeliefs.getBeliefFromId(currentBeliefId);

        // Check if the current belief is terminal
        if (terminalBeliefCallback) {
            if (auto terminal = terminalBeliefCallback(currentBelief); terminal.has_value()) {
                info.terminalBeliefValues.emplace(currentBeliefId, std::move(terminal.value()));
                continue;
            }
        }

        // Explore for each action the successors of the current belief with that action. Potentially also add rewards.
        info.exploredBeliefs.emplace(currentBeliefId, info.matrix.groups());
        auto const numActions = firstStateNextStateGenerator.getBeliefNumberOfActions(currentBelief);
        for (uint64_t localActionIndex = 0; localActionIndex < numActions; ++localActionIndex) {
            exploreNextStates(currentBelief, localActionIndex);
            info.matrix.endCurrentRow();
            if (firstStateNextStateGenerator.hasRewardModel()) {
                info.actionRewards.emplace_back(
                    storm::utility::convertNumber<BeliefMdpValueType>(firstStateNextStateGenerator.getBeliefActionReward(currentBelief, localActionIndex)));
            }
        }
        info.matrix.endCurrentRowGroup();
    }
    return true;
}

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
struct StandardDiscoverCallback {
    using InfoType = typename BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::StandardExplorationInformation;
    InfoType& info;

    StandardDiscoverCallback(InfoType& info) : info(info) {
        // Intentionally left empty
    }
    void operator()(BeliefType&& bel, typename BeliefType::ValueType&& val) {
        auto const belId = info.discoveredBeliefs.getIdOrAddBelief(std::move(bel));
        if (info.exploredBeliefs.count(belId) == 0u && info.terminalBeliefValues.count(belId) == 0u) {
            info.queue.push(belId);
        }
        info.matrix.transitions.push_back({storm::utility::convertNumber<BeliefMdpValueType>(val), belId});
    }
};

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::BeliefExploration(PomdpType const& pomdp) : firstStateNextStateGenerator(pomdp) {
    // Intentionally left empty.
}

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
typename BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::StandardExplorationInformation
BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::initializeStandardExploration(ExplorationQueueOrder const explorationQueueOrder) {
    StandardExplorationInformation info;
    info.queue.changeOrder(explorationQueueOrder);
    info.initialBeliefId = info.discoveredBeliefs.addBelief(firstStateNextStateGenerator.computeInitialBelief());
    info.queue.push(info.initialBeliefId);
    return info;
}

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
void BeliefExploration<BeliefMdpValueType, PomdpType, BeliefType>::resumeExploration(
    StandardExplorationInformation& info, TerminalBeliefCallback const& terminalBeliefCallback, TerminationCallback const& terminationCallback,
    storm::OptionalRef<std::string const> rewardModelName, storm::OptionalRef<FreudenthalTriangulationBeliefAbstraction<BeliefType>> abstraction) {
    if (rewardModelName.has_value()) {
        firstStateNextStateGenerator.setRewardModel(rewardModelName.value());
    }
    StandardDiscoverCallback<BeliefMdpValueType, PomdpType, BeliefType> discoverCallback(info);
    if (abstraction) {
        performExploration(info, firstStateNextStateGenerator.getPostAbstractionHandle(abstraction.value(), discoverCallback), terminalBeliefCallback,
                           terminationCallback);
    } else {
        performExploration(info, firstStateNextStateGenerator.getHandle(discoverCallback), terminalBeliefCallback, terminationCallback);
    }
}

template class BeliefExploration<double, storm::models::sparse::Pomdp<double>, Belief<double>>;
template class BeliefExploration<double, storm::models::sparse::Pomdp<double>, Belief<storm::RationalNumber>>;
template class BeliefExploration<storm::RationalNumber, storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<double>>;
template class BeliefExploration<storm::RationalNumber, storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<storm::RationalNumber>>;
}  // namespace storm::pomdp::beliefs
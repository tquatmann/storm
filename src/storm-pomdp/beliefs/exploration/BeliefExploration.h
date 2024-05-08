#pragma once

#include <functional>
#include <set>

#include "storm-pomdp/builder/belief-mdp/BeliefExplorationHeuristic.h"
#include "storm-pomdp/builder/belief-mdp/BeliefExplorationMatrix.h"
#include "storm-pomdp/builder/belief-mdp/BeliefExplorationMode.h"

#include "storm-pomdp/storage/beliefs/BeliefCollector.h"
#include "storm-pomdp/storage/beliefs/BeliefGenerator.h"
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"
#include "storm-pomdp/storage/beliefs/RewardBoundedBeliefSplitter.h"
#include "storm/utility/OptionalRef.h"
#include "storm/utility/vector.h"

namespace storm::pomdp::builder {

template<BeliefExplorationMode Mode, typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
class BeliefExploration {
   public:
    using BeliefId = storm::pomdp::beliefs::BeliefId;
    using StateId = storm::pomdp::beliefs::BeliefStateType;
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;

    struct ExplorationInformation {
        BeliefExplorationMatrix<BeliefMdpValueType, Mode> matrix;
        std::vector<BeliefMdpValueType> actionRewards;
        storm::pomdp::beliefs::BeliefCollector<BeliefType> collectedBeliefs;
        std::unordered_map<BeliefId, StateId> exploredBeliefs;
        BeliefId initialBelief;
    };

    struct ExpandCallback {
        ExpandCallback(ExplorationInformation& info, BeliefExplorationHeuristic<BeliefType>& explorationHeuristic)
            : info(info), explorationHeuristic(explorationHeuristic) {}
        ExplorationInformation& info;
        BeliefExplorationHeuristic<BeliefType>& explorationHeuristic;

        void operator()(BeliefType&& bel, BeliefValueType&& val) {  // TODO: add requires construct
            if constexpr (Mode == BeliefExplorationMode::Standard) {
                auto const belId = info.collectedBeliefs.getIdOrAddBelief(std::move(bel));
                if (info.exploredBeliefs.count(belId) == 0u) {
                    explorationHeuristic.discover(belId, info.collectedBeliefs.getBeliefFromId(belId));
                }
                info.matrix.transitions.push_back({storm::utility::convertNumber<BeliefMdpValueType>(val), belId});
            }
        }

        void operator()(BeliefType&& bel, std::vector<BeliefMdpValueType> const& rewardVector, BeliefValueType&& val) {  // TODO: add requires construct
            if constexpr (Mode == BeliefExplorationMode::RewardBounded) {
                auto const belId = info.collectedBeliefs.getIdOrAddBelief(std::move(bel));
                if (info.exploredBeliefs.count(belId) == 0u) {
                    explorationHeuristic.discover(belId, info.collectedBeliefs.getBeliefFromId(belId));
                }
                info.matrix.transitions.push_back({storm::utility::convertNumber<BeliefMdpValueType>(val), belId, rewardVector});
            }
        }
    };

    BeliefExploration(PomdpType const& pomdp) : beliefGenerator(pomdp), rewardBoundedBeliefSplitter(pomdp) {}

    void setRewardModelsForRewardBoundedQuery(std::vector<std::string> const& rewardModelNames) {
        rewardBoundedBeliefSplitter.setRewardModels(rewardModelNames);
    }

    void setRewardModelForObjective(std::string const& rewardModelName) {
        beliefGenerator.setRewardModel(rewardModelName);
    }

    ExplorationInformation initializeExploration(BeliefExplorationHeuristic<BeliefType>& explorationHeuristic) {
        ExplorationInformation info;
        info.initialBelief = info.collectedBeliefs.addBelief(beliefGenerator.computeInitialBelief());
        explorationHeuristic.discover(info.initialBelief, info.collectedBeliefs.getBeliefFromId(info.initialBelief));
        return info;
    }

    template<typename BeliefAbstraction>
    void exploreBelief(ExplorationInformation& info, BeliefExplorationHeuristic<BeliefType>& explorationHeuristic, BeliefId const& beliefId,
                       BeliefAbstraction const& abstraction) {
        STORM_LOG_ASSERT(info.collectedBeliefs.hasId(beliefId), "Unknown belief id");
        STORM_LOG_ASSERT(info.exploredBeliefs.count(beliefId) == 0, "Belief #" << beliefId << " already explored.");
        info.exploredBeliefs.emplace(beliefId, info.matrix.groups());
        auto belief = info.collectedBeliefs.getBeliefFromId(beliefId);  // do not take as reference since it will be invalidated when collecting more beliefs
        auto const numActions = beliefGenerator.getBeliefNumberOfActions(belief);
        ExpandCallback expandCallback(info, explorationHeuristic);
        for (uint64_t localActionIndex = 0; localActionIndex < numActions; ++localActionIndex) {
            if constexpr (Mode == BeliefExplorationMode::Standard) {
                beliefGenerator.expand(belief, localActionIndex, expandCallback, abstraction);
            } else {
                beliefGenerator.expand(rewardBoundedBeliefSplitter, belief, localActionIndex, expandCallback, abstraction);
            }
            info.matrix.endCurrentRow();
            if (beliefGenerator.hasRewardModel()) {
                info.actionRewards.emplace_back(beliefGenerator.getBeliefActionReward(belief, localActionIndex));
            }
        }
        info.matrix.endCurrentRowGroup();
    }

    template<typename BeliefAbstraction>
    bool performExploration(ExplorationInformation& info, BeliefExplorationHeuristic<BeliefType>& explorationHeuristic,
                            std::function<bool()> const& terminateExploration, BeliefAbstraction const& abstraction) {
        while (explorationHeuristic.hasNext()) {
            if (terminateExploration()) {
                return false;  // Terminate prematurely
            }
            auto currentId = explorationHeuristic.popNext();
            exploreBelief(info, explorationHeuristic, currentId, abstraction);
        }
        return true;
    }

   private:
    storm::pomdp::beliefs::BeliefGenerator<PomdpType, BeliefType> beliefGenerator;
    storm::pomdp::beliefs::RewardBoundedBeliefSplitter<BeliefMdpValueType, PomdpType, BeliefType> rewardBoundedBeliefSplitter;
};
}  // namespace storm::pomdp::builder
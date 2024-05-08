#pragma once

#include <string>
#include <vector>

#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/UnsupportedModelException.h"

namespace storm::pomdp::beliefs {

template<typename RewardValueType, typename PomdpType, typename BeliefType>
class RewardBoundedBeliefSplitter {
   public:
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;
    using RewardVectorType = std::vector<RewardValueType>;

    RewardBoundedBeliefSplitter(PomdpType const& pomdp) : pomdp(pomdp) {}

    void setRewardModel(std::string const& rewardModelName = "") {
        setRewardModels({rewardModelName});
    }

    void setRewardModels(std::vector<std::string> const& rewardModelNames) {
        actionRewardVectors.assign(pomdp.getNumberOfChoices());
        for (auto const& rewardModelName : rewardModelNames) {
            auto const& rewardModel = pomdp.getRewardModel(rewardModelName);
            STORM_LOG_THROW(!rewardModel.hasTransitionRewards(), storm::exceptions::UnsupportedModelException,
                            "POMDPs with transition rewards are currently not supported.");
            for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                for (auto const& choice : pomdp.getTransitionMatrix().getRowGroupIndices(state)) {
                    auto val = rewardModel.hasStateActionRewards() ? rewardModel.getStateActionReward(choice) : storm::utility::zero<PomdpValueType>();
                    if (rewardModel.hasStateRewards()) {
                        val += rewardModel.getStateReward(state);
                    }
                    actionRewardVectors[choice].push_back(storm::utility::convertNumber<RewardValueType>(val));
                }
            }
        }
    }

    void unsetRewardModels() {
        actionRewardVectors.clear();
    }

    std::size_t getNumberOfSetRewardModels() const {
        if (actionRewardVectors.empty()) {
            return 0;
        } else {
            return actionRewardVectors.front().size();
        }
    }

    template<typename ExpandCallback>
    void abstract(BeliefType const& belief, uint64_t localActionIndex, ExpandCallback const& callback) {
        // gather the occurring reward vectors and build the sub-beliefs
        std::unordered_map<RewardVectorType, BeliefBuilder<BeliefType>> splitBeliefs;
        belief.forEach([&localActionIndex, &splitBeliefs, this](BeliefStateType const& state, BeliefValueType const& beliefValue) {
            auto const globalActionIndex = pomdp.getTransitionMatrix().getRowGroupindices()[state] + localActionIndex;
            auto const& rewVector = actionRewardVectors[globalActionIndex];
            splitBeliefs[rewVector].addValue(state, beliefValue);
        });

        for (auto& [rewardVector, builder] : splitBeliefs) {
            builder.setObservation(belief.observation());
            if (splitBeliefs.size() == 1u) {
                // Fix the distribution to diminish numerical issues a bit
                callback(builder.build(), storm::utility::one<BeliefValueType>(), rewardVector);
            } else {
                auto val = builder.normalize();
                callback(builder.build(), std::move(val), rewardVector);
            }
        }
    }

   private:
    PomdpType const& pomdp;
    std::vector<RewardVectorType> actionRewardVectors;
};
}  // namespace storm::pomdp::beliefs
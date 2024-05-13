#pragma once

#include <string>
#include <vector>

#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/utility/constants.h"

namespace storm::pomdp::beliefs {

/*!
 * Abstracts a given belief by splitting it into sub-beliefs based on the (state-action) reward of the POMDP.
 * Multiple reward models can be set.
 * States in the support of the given belief where all action rewards are equal will be grouped together.
 * When applying this as a pre-abstraction, this intuitively means that the collected rewards are observable.
 */
template<typename RewardValueType, typename PomdpType, typename BeliefType>
class RewardBoundedBeliefSplitter {
   public:
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;
    using RewardVectorType = std::vector<RewardValueType>;

    RewardBoundedBeliefSplitter(PomdpType const& pomdp);
    void setRewardModel(std::string const& rewardModelName = "");
    void setRewardModels(std::vector<std::string> const& rewardModelNames);
    void unsetRewardModels();
    std::size_t getNumberOfSetRewardModels() const;

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
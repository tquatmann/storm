#pragma once

#include <string>
#include <vector>

#include "storm-pomdp/storage/beliefs/BeliefBuilder.h"
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {
struct NoAbstraction {};

template<typename PomdpType, typename BeliefType>
class BeliefGenerator {
   public:
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;

    BeliefGenerator(PomdpType const& pomdp) : pomdp(pomdp) {}

    void setRewardModel(std::string rewardModelName = "") {
        auto const& rewardModel = pomdp.getRewardModel(rewardModelName);
        actionRewards = rewardModel.getTotalRewardVector(pomdp.getTransitionMatrix());
    }

    bool hasRewardModel() const {
        return !actionRewards.empty();
    }

    void unsetRewardModel() {
        actionRewards.clear();
    }

    BeliefType computeInitialBelief() const {
        STORM_LOG_ASSERT(pomdp.getInitialStates().getNumberOfSetBits() == 1, "Only a single initial state is supported, but the given POMDP contains "
                                                                                 << pomdp.getInitialStates().getNumberOfSetBits() << " initial states.");
        BeliefStateType const init = *pomdp.getInitialStates().begin();
        BeliefBuilder<BeliefType> builder;
        builder.addValue(init, storm::utility::one<BeliefValueType>());
        builder.setObservation(pomdp.getObservation(init));
        return builder.build();
    }

    uint64_t getBeliefNumberOfActions(BeliefType const& belief) const {
        auto result = pomdp.getTransitionMatrix().getRowGroupSize(belief.representativeState());
        // Assert consistency with other states in the support
        STORM_LOG_ASSERT(belief.allOf([&result, this](BeliefStateType const& state, BeliefValueType const&) {
            return pomdp.getTransitionMatrix().getRowGroupSize(state) == result;
        }),
                         "Belief considers states with inconsistent number of choices.");
        return result;
    }

    PomdpValueType getBeliefActionReward(BeliefType const& belief, uint64_t const& localActionIndex) const {
        STORM_LOG_ASSERT(hasRewardModel(), "Requested a reward although no reward model was specified.");
        STORM_LOG_ASSERT(localActionIndex < getBeliefNumberOfActions(belief), "Invalid action index " << localActionIndex << ".");
        auto result = storm::utility::zero<PomdpValueType>();
        auto const& actionIndices = pomdp.getTransitionMatrix().getRowGroupIndices();
        belief.forEach([&localActionIndex, &actionIndices, &result, this](BeliefStateType const& state, BeliefValueType const& val) {
            uint64_t const actionIndex = actionIndices[state] + localActionIndex;
            result += storm::utility::convertNumber<PomdpValueType>(val) * this->actionRewards[actionIndex];
        });
        return result;
    }

    template<typename ExpandCallback>
    void expand(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& probabilityFactor, ExpandCallback& callback) {
        expand(belief, localActionIndex, probabilityFactor, callback, NoAbstraction{});
    }

    template<typename ExpandCallback, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, ExpandCallback& callback, PostAbstraction& abstraction) {
        expand(belief, localActionIndex, storm::utility::one<BeliefValueType>(), callback, abstraction);
    }

    template<typename ExpandCallback, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& probabilityFactor, ExpandCallback& callback,
                PostAbstraction& abstraction) {
        // Find the probability we go to each observation
        std::unordered_map<BeliefObservationType, BeliefValueType> successorObservations;
        belief.forEach([&localActionIndex, &successorObservations, this](BeliefStateType const& state, BeliefValueType const& beliefValue) {
            for (auto const& pomdpTransition : pomdp.getTransitionMatrix().getRow(state, localActionIndex)) {
                if (!storm::utility::isZero(pomdpTransition.getValue())) {
                    auto const obs = pomdp.getObservation(pomdpTransition.getColumn());
                    auto const val = beliefValue * storm::utility::convertNumber<BeliefValueType>(pomdpTransition.getValue());
                    if (auto [insertionIt, inserted] = successorObservations.emplace(obs, val); !inserted) {
                        insertionIt->second += val;
                    }
                }
            }
        });

        // Adjust the distribution to diminish numerical inacuracies a bit
        if constexpr (!storm::NumberTraits<BeliefValueType>::IsExact || !storm::NumberTraits<PomdpValueType>::IsExact) {
            if (successorObservations.size() == 1) {
                successorObservations.begin()->second = storm::utility::one<BeliefValueType>();
            }
        }

        // For each successor observation we build the successor belief
        for (auto const& successorObsValue : successorObservations) {
            BeliefBuilder<BeliefType> builder;
            builder.setObservation(successorObsValue.first);
            belief.forEach([&builder, &localActionIndex, &successorObsValue, this](BeliefStateType const& state, BeliefValueType const& beliefValue) {
                for (auto const& pomdpTransition : pomdp.getTransitionMatrix().getRow(state, localActionIndex)) {
                    if (pomdp.getObservation(pomdpTransition.getColumn()) == successorObsValue.first) {
                        BeliefValueType const prob =
                            beliefValue * storm::utility::convertNumber<BeliefValueType>(pomdpTransition.getValue()) / successorObsValue.second;
                        builder.addValue(pomdpTransition.getColumn(), prob);
                    }
                }
            });
            if constexpr (std::is_same_v<std::remove_cv_t<PostAbstraction>, NoAbstraction>) {
                callback(static_cast<BeliefValueType>(successorObsValue.second * probabilityFactor), builder.build());
            } else {
                abstraction.abstract(static_cast<BeliefValueType>(successorObsValue.second * probabilityFactor), builder.build(), callback);
            }
        }
    }

    template<typename ExpandCallback, typename PreAbstraction, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, ExpandCallback& callback, PreAbstraction& preAbstraction,
                PostAbstraction& postAbstraction) {
        if constexpr (PreAbstraction::CallbackHasData) {
            using CbData = typename PreAbstraction::CallbackDataType;
            preAbstraction.abstract(
                belief, localActionIndex,
                [this, &localActionIndex, &callback, &postAbstraction](CbData const& data, BeliefValueType const& preVal, BeliefType const& preBel) {
                    expand(preBel, localActionIndex, preVal, [&data, &callback](auto&&... args) { callback(data, args...); }, postAbstraction);
                });
        } else {
            preAbstraction.abstract(belief, localActionIndex,
                                    [this, &localActionIndex, &callback, &postAbstraction](BeliefValueType const& preVal, BeliefType&& preBel) {
                                        expand(preBel, localActionIndex, preVal, callback, postAbstraction);
                                    });
        }
    }

   private:
    PomdpType const& pomdp;
    std::vector<PomdpValueType> actionRewards;
};
}  // namespace storm::pomdp::beliefs
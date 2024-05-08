#pragma once

#include <string>
#include <vector>

#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {
struct NoAbstractionStruct {};
using NoAbstraction = NoAbstractionStruct const;
inline constexpr NoAbstraction noAbstraction;

template<typename PomdpType, typename BeliefType>
class FirstStateNextStateGenerator {
   public:
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;

    FirstStateNextStateGenerator(PomdpType const& pomdp) : pomdp(pomdp) {}

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

    template<typename ExpandCallback, typename PreAbstraction = NoAbstraction, typename PostAbstraction = NoAbstraction>
    struct Expander {
        PomdpType const& pomdp;
        ExpandCallback& callback;
        PreAbstraction& preAbstraction;
        PostAbstraction& postAbstraction;

        void operator()(BeliefType const& belief, uint64_t localActionIndex,
                        BeliefValueType const& transitionProbability = storm::utility::one<BeliefValueType>()) {
            applyPreAbstraction(belief, localActionIndex, transitionProbability);
        }

       private:
        template<typename... CallBackArgs>
        void applyPreAbstraction(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& transitionProbability,
                                 CallBackArgs const&... additionalCallbackArgs) {
            if constexpr (std::is_same_v<std::remove_cv_t<PreAbstraction>, std::remove_cv_t<NoAbstraction>>) {
                computeSuccessorBeliefs(belief, localActionIndex, transitionProbability, additionalCallbackArgs...);
            } else {
                preAbstraction.abstract(
                    belief, localActionIndex, transitionProbability,
                    [this, &localActionIndex, &additionalCallbackArgs...](BeliefType&& preBel, BeliefValueType&& preVal, auto const&... extraData) {
                        computeSuccessorBeliefs(preBel, localActionIndex, preVal, additionalCallbackArgs..., extraData...);
                    });
            }
            //        } else if constexpr (PreAbstraction::CallbackHasData) {
            //
            //            using CbData = typename PreAbstraction::CallbackDataType;
            //            preAbstraction.abstract(
            //                belief, localActionIndex, transitionProbability,
            //                [this, &localActionIndex, &additionalCallbackArgs...](CbData const& data, BeliefValueType const& preVal, BeliefType const& preBel)
            //                {
            //                    computeSuccessorBeliefs(preBel, localActionIndex, preVal, data, additionalCallbackArgs...);
            //                });
            //        }
        }

        template<typename... CallBackArgs>
        void computeSuccessorBeliefs(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& transitionProbability,
                                     CallBackArgs const&... additionalCallbackArgs) {
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
                applyPostAbstraction(builder.build(), static_cast<BeliefValueType>(successorObsValue.second * transitionProbability),
                                     additionalCallbackArgs...);
            }
        }

        template<typename... CallBackArgs>
        void applyPostAbstraction(BeliefType&& belief, BeliefValueType&& transitionProbability, CallBackArgs const&... additionalCallbackArgs) {
            if constexpr (std::is_same_v<std::remove_cv_t<PostAbstraction>, std::remove_cv_t<NoAbstraction>>) {
                callback(std::move(belief), std::move(transitionProbability), additionalCallbackArgs...);
            } else {
                postAbstraction.abstract(
                    std::move(belief), std::move(transitionProbability),
                    [this, &additionalCallbackArgs...](BeliefType&& postBel, BeliefValueType&& postVal, auto const&... postAbstractionData) {
                        callback(std::move(postBel), std::move(postVal), additionalCallbackArgs..., postAbstractionData...);
                    });
            }
        }
    };

    template<typename ExpandCallback>
    void expand(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& probabilityFactor, ExpandCallback& callback) {
        Expander<ExpandCallback, NoAbstraction, NoAbstraction> expander{this->pomdp, callback, noAbstraction, noAbstraction};
        expander(belief, localActionIndex, probabilityFactor);
    }

    template<typename ExpandCallback, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, ExpandCallback& callback, PostAbstraction& postAbstraction) {
        Expander<ExpandCallback, NoAbstraction, PostAbstraction> expander{this->pomdp, callback, noAbstraction, postAbstraction};
        expander(belief, localActionIndex);
    }

    template<typename ExpandCallback, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, BeliefValueType const& probabilityFactor, ExpandCallback& callback,
                PostAbstraction& postAbstraction) {
        Expander<ExpandCallback, NoAbstraction, PostAbstraction> expander{this->pomdp, callback, noAbstraction, postAbstraction};
        expander(belief, localActionIndex, probabilityFactor);
        /* old
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
         */
    }

    template<typename ExpandCallback, typename PreAbstraction, typename PostAbstraction>
    void expand(BeliefType const& belief, uint64_t localActionIndex, ExpandCallback& callback, PreAbstraction& preAbstraction,
                PostAbstraction& postAbstraction) {
        Expander<ExpandCallback, PreAbstraction, PostAbstraction> expander{this->pomdp, callback, preAbstraction, postAbstraction};
        expander(belief, localActionIndex);
        /* old
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
         */
    }

   private:
    PomdpType const& pomdp;
    std::vector<PomdpValueType> actionRewards;
};
}  // namespace storm::pomdp::beliefs
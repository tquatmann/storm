#pragma once
#include "storm-pomdp/beliefs/utility/types.h"

#include "storm-pomdp/beliefs/utility/BeliefNumerics.h"

#include "storm/utility/constants.h"

#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
class BeliefBuilder {
    using ValueType = typename BeliefType::ValueType;

   public:
    void reserve(uint64_t size) {
        data.reserve(size);
    }

    void addValue(BeliefStateType const& state, ValueType const& value) {
        STORM_LOG_ASSERT(value > storm::utility::zero<ValueType>() && value <= storm::utility::one<ValueType>(), "Invalid belief value " << value << ".");
        if (auto [insertionIt, success] = data.emplace(state, value); !success) {
            insertionIt->second += value;
        }
    }

    void setObservation(BeliefObservationType const& observation) {
        obs = observation;
    }

    ValueType normalize() {
        auto sum = storm::utility::zero<ValueType>();
        for (auto const& [state, value] : data) {
            STORM_LOG_ASSERT(value >= storm::utility::zero<ValueType>(), "Invalid value: " << value);
            sum += value;
        }
        if (!storm::utility::isOne(sum)) {
            for (auto& [state, value] : data) {
                value /= sum;
            }
        }
        return sum;
    }

    BeliefType build() {
        data.shrink_to_fit();
        if constexpr (!storm::NumberTraits<ValueType>::IsExact) {
            if (data.size() == 1 && BeliefNumerics<ValueType>::isOne(data.begin()->second)) {
                // If the distribution consists of only one entry and its value is sufficiently close to 1, make it exactly 1 to avoid numerical problems
                data.begin()->second = storm::utility::one<ValueType>();
            }
        }
        STORM_LOG_ASSERT(assertBelief(), "Trying to build invalid belief.");
        return BeliefType{std::move(data), std::move(obs)};
    }

   private:
    bool assertBelief() {
        using ValueType = typename BeliefType::ValueType;
        if (obs == InvalidObservation) {
            STORM_LOG_ERROR("Observation of Belief not set.");
            return false;
        }
        auto sum = storm::utility::zero<ValueType>();
        for (auto const& [state, value] : data) {
            if (value <= storm::utility::zero<ValueType>()) {
                STORM_LOG_ERROR("Non-positive belief value " << value << " at state " << state << ".");
                return false;
            }
            if (!BeliefNumerics<ValueType>::lessOrEqual(value, storm::utility::one<ValueType>())) {
                STORM_LOG_ERROR("Invalid belief value " << value << " at state " << state << ".");
                return false;
            }
            sum += value;
        }
        if (!BeliefNumerics<ValueType>::lessOrEqual(sum, storm::utility::one<ValueType>())) {
            STORM_LOG_ERROR("belief value sum " << sum << " is larger than 1 (sum-1=" << sum - storm::utility::one<ValueType>() << ").");
            return false;
        }
        if (!BeliefNumerics<ValueType>::lessOrEqual(storm::utility::one<ValueType>(), sum)) {
            STORM_LOG_ERROR("belief value sum " << sum << " is smaller than 1 (1-sum=" << storm::utility::one<ValueType>() - sum << ").");
            return false;
        }
        return true;
    }

    BeliefFlatMap<ValueType> data;
    BeliefObservationType obs{InvalidObservation};
};

}  // namespace storm::pomdp::beliefs
#pragma once

#include <sstream>

#include <boost/functional/hash.hpp>

#include "storm/utility/NumberTraits.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/utility/BeliefNumerics.h"
#include "storm-pomdp/beliefs/utility/types.h"

namespace storm::pomdp::beliefs {

template<typename ValueTypeArg>
class Belief {
   public:
    using ValueType = ValueTypeArg;
    friend class BeliefBuilder<Belief<ValueType>>;

    Belief() = delete;
    Belief(Belief const& other) = default;
    Belief(Belief&& other) = default;
    Belief<ValueType>& operator=(Belief const& other) = default;
    Belief<ValueType>& operator=(Belief&& other) = default;

    template<typename FunctionType>
    void forEach(FunctionType const& f) const {
        for (auto const& [state, value] : data) {
            f(state, value);
        }
    }

    template<typename FunctionType>
    bool allOf(FunctionType const& f) const {
        for (auto const& [state, value] : data) {
            if (!f(state, value)) {
                return false;
            }
        }
        return true;
    }

    template<typename FunctionType>
    void forEachStateInSupport(FunctionType const& f) const {
        for (auto const& [state, value] : data) {
            f(state);
        }
    }
    template<typename FunctionType>
    void forEachCombine(Belief<ValueType> const& other, FunctionType const& f, bool considerOnlyThisSupport = false) const {
        static_assert(BeliefFlatMapIsOrdered);
        auto const zero = storm::utility::zero<ValueType>();
        auto it2 = other.data.cbegin();
        auto const it2End = data.cend();
        if (considerOnlyThisSupport) {
            for (auto const& [state1, value1] : data) {
                while (it2 != it2End && it2->first < state1) {
                    ++it2;
                }
                if (it2 == it2End || it2->first > state1) {
                    f(state1, value1, zero);
                } else {
                    STORM_LOG_ASSERT(state1 == it2->first, "Unexpected state.");
                    f(state1, value1, it2->second);
                }
            }
        } else {
            auto it1 = data.cbegin();
            auto const it1End = data.cend();
            while (it1 != it1End && it2 != it2End) {
                auto const state1 = it1->first;
                auto const state2 = it2->first;
                if (state1 == state2) {
                    f(state1, it1->second, it2->second);
                    ++it1;
                    ++it2;
                } else if (state1 < state2) {
                    f(state1, it1->second, zero);
                    ++it1;
                } else {
                    STORM_LOG_ASSERT(state2 < state1, "unexpected states.");
                    f(state2, zero, it2->second);
                    ++it2;
                }
            }
            for (; it1 != it1End; ++it1) {
                f(it1->first, it1->second, zero);
            }
            for (; it2 != it2End; ++it2) {
                f(it2->first, zero, it2->second);
            }
        }
    }

    std::size_t size() const {
        return data.size();
    }

    BeliefStateType representativeState() const {
        STORM_LOG_ASSERT(!data.empty(), "Empty belief");
        return data.begin()->first;
    }

    BeliefObservationType const& observation() const {
        return obs;
    }

    bool operator==(Belief const& other) const {
        if (obs != other.obs) {
            return false;
        }
        if (data.size() != other.size()) {
            return false;
        }
        static_assert(BeliefFlatMapIsOrdered);
        auto secondIt = other.data.cbegin();
        for (auto const& [state, value] : data) {
            if (state != secondIt->first) {
                return false;
            }
            if (!BeliefNumerics<ValueType>::equal(value, secondIt->second)) {
                return false;
            }
            ++secondIt;
        }
        return true;
    }

    std::string toString(bool convertToDouble = true) const {
        std::stringstream ss;
        ss << "Belief{ obs:" << obs;
        if (convertToDouble) {
            forEach([&ss](auto const& state, auto const& val) { ss << ", " << state << ":" << storm::utility::convertNumber<double>(val); });
        } else {
            forEach([&ss](auto const& state, auto const& val) { ss << ", " << state << ":" << val; });
        }
        ss << " }";
        return ss.str();
    }

    template<typename SummandsType>
    SummandsType getWeightedSum(std::vector<SummandsType> const& summands) {
        auto sum = storm::utility::zero<SummandsType>();
        forEach([&sum, &summands](auto const& state, auto const& val) {
            STORM_LOG_ASSERT(state < summands.size(), "State " << state << " is out of range for the given summands.");
            sum += summands[state] * val;
        });
        return sum;
    }

    struct BeliefHash {
        std::size_t operator()(Belief const& belief) const {
            auto seed = static_cast<std::size_t>(belief.obs);
            if constexpr (storm::NumberTraits<ValueType>::IsExact) {
                boost::hash_combine(seed, belief.data);
            } else {
                static_assert(BeliefFlatMapIsOrdered);
                belief.forEach([&seed](auto const& state, auto const& val) {
                    boost::hash_combine(seed, state);
                    boost::hash_combine(seed, BeliefNumerics<ValueType>::valueForHash(val));
                });
            }
            return seed;
        }
    };

   private:
    Belief(BeliefFlatMap<ValueType>&& data, BeliefObservationType&& obs) : data(std::move(data)), obs(std::move(obs)) {}
    BeliefFlatMap<ValueType> const data;
    BeliefObservationType const obs;
};

}  // namespace storm::pomdp::beliefs
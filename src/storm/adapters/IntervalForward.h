#pragma once

#include "storm/adapters/RationalNumberForward.h"

namespace carl {
template<typename Number>
class Interval;
}

namespace storm {

/*!
 * Interval type
 */
typedef carl::Interval<double> Interval;
typedef carl::Interval<storm::RationalNumber> RationalInterval;

namespace detail {
template<typename ValueType>
struct IntervalMetaProgrammingHelper {
    using BaseType = ValueType;
    // no interval type in the generic setting; fall back to the type itself.
    using IntervalType = ValueType;
    static const bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<double> {
    using BaseType = double;
    using IntervalType = Interval;
    static const bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<storm::RationalNumber> {
    using BaseType = storm::RationalNumber;
    using IntervalType = RationalInterval;
    static const bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<Interval> {
    using BaseType = double;
    using IntervalType = Interval;
    static const bool isInterval = true;
};
template<>
struct IntervalMetaProgrammingHelper<RationalInterval> {
    using BaseType = storm::RationalNumber;
    using IntervalType = RationalInterval;
    static const bool isInterval = true;
};
}  // namespace detail

/*!
 * Helper to check if a type is an interval
 */
template<typename ValueType>
constexpr bool IsIntervalType = detail::IntervalMetaProgrammingHelper<ValueType>::isInterval;

/*!
 * Helper to access the type in which interval boundaries are stored.
 * Yields the type identity if the given type is not an interval
 */
template<typename ValueType>
using IntervalBaseType = typename detail::IntervalMetaProgrammingHelper<ValueType>::BaseType;

/*!
 * Helper to access the associated interval type of the given ValueType.
 * Yields the type identity if the given type is already an interval
 */
template<typename ValueType>
using IntervalType = typename detail::IntervalMetaProgrammingHelper<ValueType>::IntervalType;

}  // namespace storm

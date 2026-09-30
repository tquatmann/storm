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
    // Interval no interval type available for the default helper (e.g. used for RationalFunction)
    static constexpr bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<double> {
    using BaseType = double;
    using IntervalType = Interval;
    static constexpr bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<storm::RationalNumber> {
    using BaseType = storm::RationalNumber;
    using IntervalType = RationalInterval;
    static constexpr bool isInterval = false;
};
template<>
struct IntervalMetaProgrammingHelper<Interval> {
    using BaseType = double;
    using IntervalType = Interval;
    static constexpr bool isInterval = true;
};
template<>
struct IntervalMetaProgrammingHelper<RationalInterval> {
    using BaseType = storm::RationalNumber;
    using IntervalType = RationalInterval;
    static constexpr bool isInterval = true;
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
 * Helper to access the interval type whose bounds are of the given type, e.g., storm::Interval for double.
 * Yields the type identity if the given type already is an interval type and is not defined if there is no such interval type.
 */
template<typename ValueType>
using IntervalType = typename detail::IntervalMetaProgrammingHelper<ValueType>::IntervalType;
}  // namespace storm

#include "storm-pomdp/storage/beliefs/BeliefNumerics.h"

#include <cmath>
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/utility/Constants.h"
#include "storm/utility/NumberTraits.h"

namespace storm::pomdp::beliefs {

namespace detail {
/*!
 * Used to decide whether two belief values are equal.
 * The interpretation is that two values x,y are equal if round(x*2^NumericPrecision)==round(y*2^NumericPrecision).
 * Thus, a high value means that we are more precise.
 * This is only relevant if beliefs are represented using an inexact data type like double
 */
constexpr double NumericPrecisionFactor = 1e15;

double reprValue(double const& val) {
    return round(val * NumericPrecisionFactor);
}

storm::RationalNumber reprValue(storm::RationalNumber const& val) {
    return val;
}
}  // namespace detail

template<typename ValueType>
bool BeliefNumerics<ValueType>::lessOrEqual(ValueType const& lhs, ValueType const& rhs) {
    return detail::reprValue(lhs) <= detail::reprValue(rhs);
}

template<typename ValueType>
bool BeliefNumerics<ValueType>::equal(ValueType const& lhs, ValueType const& rhs) {
    return detail::reprValue(lhs) == detail::reprValue(rhs);
}

template<typename ValueType>
bool BeliefNumerics<ValueType>::isZero(const ValueType& val) {
    return storm::utility::isZero(detail::reprValue(val));
}

template<typename ValueType>
bool BeliefNumerics<ValueType>::isOne(const ValueType& val) {
    return storm::utility::isOne(detail::reprValue(val));
}

template<typename ValueType>
ValueType BeliefNumerics<ValueType>::valueForHash(const ValueType& val) {
    return detail::reprValue(val);
}

template class BeliefNumerics<double>;
template class BeliefNumerics<storm::RationalNumber>;

}  // namespace storm::pomdp::beliefs
#pragma once

#include "storm/adapters/IntervalForward.h"
#include "storm/exceptions/InvalidOperationException.h"
#include "storm/storage/Distribution.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"

namespace storm {
namespace utility {
namespace interval {

/*!
 * Takes the given interval and computes the intersection with [0, 1] with respect to the precision of the given comparator.
 * @param interval
 * @throws InvalidOperationException if the intersection is empty.
 * @return probabilistic interval I s.t. I = interval \\cap [0, 1]
 */
template<typename ValueType>
ValueType makeIntervalProbabilistic(ValueType interval, ConstantsComparator<IntervalBaseType<ValueType>> const& comparator) {
    using BaseType = IntervalBaseType<ValueType>;

    auto const zero = storm::utility::zero<BaseType>();
    auto const one = storm::utility::one<BaseType>();

    auto lowerBound = interval.lower();
    auto upperBound = interval.upper();

    // Correct tiny numerical drift at interval bounds.
    if (comparator.isEqual(lowerBound, zero)) {
        lowerBound = zero;
    }
    if (comparator.isEqual(upperBound, zero)) {
        upperBound = zero;
    }
    if (comparator.isEqual(lowerBound, one)) {
        lowerBound = one;
    }
    if (comparator.isEqual(upperBound, one)) {
        upperBound = one;
    }

    // Compute I \cap [0, 1].
    auto intersectedLower = std::max(lowerBound, zero);
    auto intersectedUpper = std::min(upperBound, one);

    STORM_LOG_THROW(!comparator.isLess(intersectedUpper, intersectedLower), storm::exceptions::InvalidOperationException,
                    "Interval " << interval << " has empty intersection with [0, 1]. "
                                << "After tolerance-based bounds correction: [" << lowerBound << ", " << upperBound << "], intersection bounds: ["
                                << intersectedLower << ", " << intersectedUpper << "].");

    return ValueType(intersectedLower, carl::BoundType::WEAK, intersectedUpper, carl::BoundType::WEAK);
}

/*!
 * Takes the given interval and computes the intersection with [0, 1].
 * @param interval
 * @throws InvalidOperationException if the intersection is empty.
 * @return probabilistic interval I s.t. I = interval \\cap [0, 1]
 */
template<typename ValueType>
ValueType makeIntervalProbabilistic(ValueType interval) {
    return makeIntervalProbabilistic<ValueType>(interval,
                                                ConstantsComparator<IntervalBaseType<ValueType>>(storm::utility::zero<IntervalBaseType<ValueType>>()));
}

/*!
 * Returns the given distribution where each interval entry I gets turned into a probabilistic interval I', i.e., I' = I \cap [0, 1], with respect to the
 * precision of the given comparator.
 * @param distribution
 * @throws InvalidOperationException if any of the intervals contained in the distribution are empty after being made probabilistic.
 * @return distribution with probabilistic interval entries
 */
template<typename ValueType>
storage::Distribution<ValueType, storm::storage::sparse::state_type> makeDistributionProbabilistic(
    storage::Distribution<ValueType, storm::storage::sparse::state_type> distribution, ConstantsComparator<IntervalBaseType<ValueType>> const& comparator) {
    auto intersectedDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
    intersectedDistribution.reserve(distribution.size());

    for (auto it = distribution.begin(); it != distribution.end(); ++it) {
        auto key = it->first;
        auto interval = it->second;

        intersectedDistribution.addProbability(key, storm::utility::interval::makeIntervalProbabilistic(interval, comparator));
    }

    return intersectedDistribution;
}

/*!
 * Returns the given distribution where each interval entry I gets turned into a probabilistic interval I', i.e., I' = I \cap [0, 1].
 * @param distribution
 * @throws InvalidOperationException if any of the intervals contained in the distribution are empty after being made probabilistic.
 * @return distribution with probabilistic interval entries
 */
template<typename ValueType>
storage::Distribution<ValueType, storm::storage::sparse::state_type> makeDistributionProbabilistic(
    storage::Distribution<ValueType, storm::storage::sparse::state_type> distribution) {
    return makeDistributionProbabilistic(distribution, ConstantsComparator<ValueType>(storm::utility::zero<IntervalBaseType<ValueType>>()));
}

}  // namespace interval
}  // namespace utility
}  // namespace storm
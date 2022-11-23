#pragma once

#include <optional>

#include "storm/modelchecker/prctl/helper/BaierUpperRewardBoundsComputer.h"
#include "storm/modelchecker/prctl/helper/DsMpiUpperRewardBoundsComputer.h"
#include "storm/solver/AbstractEquationSolver.h"
#include "storm/solver/OptimizationDirection.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

#include "storm/exceptions/InvalidOperationException.h"

namespace storm::solver::helper {

template<typename ValueType>
void computeLowerUpperBounds(AbstractEquationSolver<ValueType>& solver, storm::storage::SparseMatrix<ValueType> const& matrix,
                             std::vector<ValueType> const& offsets, std::vector<ValueType> const& exitProbabilities, bool reqLower, bool reqUpper,
                             bool stronglyConnected = false, std::optional<storm::OptimizationDirection> dir = {}) {
    if constexpr (std::is_same_v<ValueType, storm::RationalFunction>) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Operation not allowed.");
    } else {
        STORM_LOG_ASSERT(matrix.hasTrivialRowGrouping() || dir.has_value(), "Optimization direction needs to be given if matrix is nondeterministic.");
        storm::utility::Extremum<storm::OptimizationDirection::Minimize, ValueType> lowerBound;
        storm::utility::Extremum<storm::OptimizationDirection::Maximize, ValueType> upperBound;

        // Try to compute a lower and an upper bound the "easy" way
        bool canLower = true;
        bool canUpper = true;
        for (uint64_t row = 0; row < offsets.size(); ++row) {
            if (storm::utility::isZero(exitProbabilities[row])) {
                if (canLower && offsets[row] < storm::utility::zero<ValueType>()) {
                    // We cannot compute a lower bound
                    canLower = false;
                }
                if (canUpper && offsets[row] > storm::utility::zero<ValueType>()) {
                    // We cannot compute an upper bound
                    canUpper = false;
                }
                if (!canUpper && !canLower) {
                    break;
                }
            } else {
                ValueType rowValue = offsets[row] / exitProbabilities[row];
                if (canUpper) {
                    upperBound &= rowValue;
                }
                if (canLower) {
                    lowerBound &= rowValue;
                }
            }
        }

        if (canUpper) {
            solver.setUpperBound(*upperBound);
        }
        if (canLower) {
            solver.setLowerBound(*lowerBound);
        }

        // We might have to invoke the respective reward bound computers.
        std::vector<ValueType> tmpOffsets;
        if (reqUpper && !canUpper) {
            bool hasNegativeValues = !canLower || *lowerBound < storm::utility::zero<ValueType>();
            if (hasNegativeValues) {
                tmpOffsets.resize(offsets.size());
                storm::utility::vector::applyPointwise(offsets, tmpOffsets, [](ValueType const& v) { return std::max(storm::utility::zero<ValueType>(), v); });
            }
            if (dir.has_value() && maximize(*dir)) {
                std::function<uint64_t(uint64_t)> stateToScc;
                if (stronglyConnected) {
                    stateToScc = [](uint64_t) -> uint64_t { return 0; };
                }
                solver.setUpperBound(storm::modelchecker::helper::BaierUpperRewardBoundsComputer<ValueType>(matrix, hasNegativeValues ? tmpOffsets : offsets,
                                                                                                            exitProbabilities, stateToScc)
                                         .computeUpperBound());
            } else {
                solver.setUpperBounds(storm::modelchecker::helper::DsMpiMdpUpperRewardBoundsComputer<ValueType>(
                                          matrix, hasNegativeValues ? tmpOffsets : offsets, exitProbabilities)
                                          .computeUpperBounds());
            }
        }
        if (reqLower && !canLower) {
            // For lower bounds we actually compute upper bounds for the negated rewards.
            // We therefore need tmpOffsets in anyway.
            tmpOffsets.resize(offsets.size());
            storm::utility::vector::applyPointwise(offsets, tmpOffsets,
                                                   [](ValueType const& v) { return std::max<ValueType>(storm::utility::zero<ValueType>(), -v); });
            if (dir.has_value() && minimize(*dir)) {
                std::function<uint64_t(uint64_t)> stateToScc;
                if (stronglyConnected) {
                    stateToScc = [](uint64_t) -> uint64_t { return 0; };
                }
                solver.setLowerBound(-storm::modelchecker::helper::BaierUpperRewardBoundsComputer<ValueType>(matrix, tmpOffsets, exitProbabilities, stateToScc)
                                          .computeUpperBound());
            } else {
                auto lowerBounds =
                    storm::modelchecker::helper::DsMpiMdpUpperRewardBoundsComputer<ValueType>(matrix, tmpOffsets, exitProbabilities).computeUpperBounds();
                storm::utility::vector::applyPointwise(lowerBounds, lowerBounds, [](ValueType const& v) { return -v; });
                solver.setLowerBounds(std::move(lowerBounds));
            }
        }
    }
}
}  // namespace storm::solver::helper
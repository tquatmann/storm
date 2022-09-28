
#pragma once
#include <functional>
#include <optional>
#include <vector>

#include "storm/solver/SolverStatus.h"
#include "storm/solver/helper/ValueIterationHelper.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/KwekMehlhorn.h"
#include "storm/utility/vector.h"

namespace storm {
class Environment;

namespace solver::helper {

template<typename TargetValueType, typename ExactValueType, typename ImpreciseValueType, bool TrivialRowGrouping = false>
class RationalSearchHelper {
   public:
    static const bool IsTargetExact = std::is_same_v<TargetValueType, ExactValueType>;

    RationalSearchHelper(std::shared_ptr<ValueIterationOperator<ExactValueType, TrivialRowGrouping>> exactViOperator,
                         std::shared_ptr<ValueIterationOperator<ImpreciseValueType, TrivialRowGrouping>> impreciseViOperator)
        : _exactOperator(exactViOperator), _impreciseOperator(impreciseViOperator) {
        // Intentionally left empty
    }

    enum class RSResult { InProgress, Converged, PrecisionExceeded };

    template<typename ValueType, storm::OptimizationDirection Dir>
    RSResult sharpen(uint64_t precision, std::vector<ValueType> const& operand, std::vector<ExactValueType> const& exactOffsets,
                     std::vector<TargetValueType>& target) {
        RSBackend<Dir> backend;
        auto& sharpOperand = _exactOperator->allocateAuxiliaryVector(operand.size());

        for (uint64_t p = 0; p <= precision; ++p) {
            // If we need an exact computation but are currently using imprecise values, we might need to consider switching to precise values
            if (std::is_same_v<ExactValueType, TargetValueType> && std::is_same_v<ValueType, ImpreciseValueType> &&
                p > std::numeric_limits<ImpreciseValueType>::max_digits10) {
                _exactOperator->freeAuxiliaryVector();
                return RSResult::PrecisionExceeded;
            }
            storm::utility::kwek_mehlhorn::sharpen(p, operand, sharpOperand);

            if (_exactOperator->applyInPlace(sharpOperand, exactOffsets, backend)) {
                // Put the solution into the target vector
                if constexpr (std::is_same_v<ExactValueType, TargetValueType>) {
                    target.swap(sharpOperand);
                } else {
                    target.resize(sharpOperand.size());
                    storm::utility::vector::convertNumericVector(sharpOperand, target);
                }
                _exactOperator->freeAuxiliaryVector();
                return RSResult::Converged;
            }
        }
        _exactOperator->freeAuxiliaryVector();
        return RSResult::InProgress;
    }

    template<storm::OptimizationDirection Dir>
    class RSBackend {
       public:
        void startNewIteration() {
            _allEqual = true;
        }

        void firstRow(ExactValueType&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            _best = std::move(value);
        }

        void nextRow(ExactValueType&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            _best &= value;
        }

        void applyUpdate(ExactValueType& currValue, [[maybe_unused]] uint64_t rowGroup) {
            if (currValue != *_best) {
                _allEqual = false;
            }
        }

        void endOfIteration() const {
            // intentionally left empty.
        }

        bool converged() const {
            return _allEqual;
        }

        bool constexpr abort() const {
            return !_allEqual;
        }

       private:
        storm::utility::Extremum<Dir, ExactValueType> _best;
        bool _allEqual{true};
    };

    template<typename ValueType>
    auto& getOperator() {
        if constexpr (std::is_same_v<ValueType, ExactValueType>) {
            return _exactOperator;
        } else {
            return _impreciseOperator;
        }
    }

    template<typename ValueType, storm::OptimizationDirection Dir>
    auto RS(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, ExactValueType precision,
            std::vector<ExactValueType> const& exactOffsets, std::vector<TargetValueType>& target,
            std::function<SolverStatus(SolverStatus const&)> const& iterationCallback) {
        ValueIterationHelper<ValueType, TrivialRowGrouping> viHelper(getOperator<ValueType>());
        SolverStatus status{SolverStatus::InProgress};
        RSResult result{RSResult::InProgress};
        while (status == SolverStatus::InProgress) {
            auto viStatus =
                viHelper.template VI<Dir, false>(operand, offsets, numIterations, storm::utility::convertNumber<ValueType>(precision), iterationCallback);

            // Compute maximal precision until which to sharpen.
            auto p = storm::utility::convertNumber<uint64_t>(
                storm::utility::ceil(storm::utility::log10<ExactValueType>(storm::utility::one<ExactValueType>() / precision)));

            // check if the sharpened vector is the desired solution.
            result = sharpen<ValueType, Dir>(p, operand, exactOffsets, target);
            switch (result) {
                case RSResult::Converged:
                    status = SolverStatus::Converged;
                    break;
                case RSResult::InProgress:
                    if (viStatus != SolverStatus::Converged) {
                        status = viStatus;
                    } else {
                        // Increase the precision.
                        precision /= storm::utility::convertNumber<ExactValueType>(static_cast<uint64_t>(10));
                    }
                    break;
                case RSResult::PrecisionExceeded:
                    status = SolverStatus::Aborted;
                    break;
            }
        }
        return std::pair(result, status);
    }

    template<storm::OptimizationDirection Dir>
    auto RS(std::vector<TargetValueType>& operand, std::vector<TargetValueType> const& offsets, uint64_t& numIterations, TargetValueType const& precision,
            std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        if constexpr (std::is_same_v<TargetValueType, ExactValueType>) {
            // We need an exact solution
            // We first try to solve the problem using imprecise values and fall back to exact values if needed.
            auto& impreciseOperand = _impreciseOperator->allocateAuxiliaryVector(operand.size());
            storm::utility::vector::convertNumericVector(operand, impreciseOperand);
            auto impreciseOffsets = storm::utility::vector::convertNumericVector<ImpreciseValueType>(offsets);
            auto [result, status] =
                RS<ImpreciseValueType, Dir>(impreciseOperand, impreciseOffsets, numIterations, precision, offsets, operand, iterationCallback);
            if (result != RSResult::PrecisionExceeded) {
                return status;
            }
            STORM_LOG_WARN("Precision of value type was exceeded, trying to recover by switching to rational arithmetic.");
            storm::utility::vector::convertNumericVector(impreciseOperand, operand);
            return RS<TargetValueType, Dir>(operand, offsets, numIterations, precision, offsets, operand, iterationCallback).second;
        } else {
            // We only try with the inexact type
            auto exactOffsets = storm::utility::vector::convertNumericVector<ExactValueType>(offsets);
            return RS<TargetValueType, Dir>(operand, offsets, numIterations, storm::utility::convertNumber<ExactValueType>(precision), exactOffsets, operand,
                                            iterationCallback)
                .second;
        }
    }

    auto RS(std::vector<TargetValueType>& operand, std::vector<TargetValueType> const& offsets, uint64_t& numIterations, TargetValueType const& precision,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        STORM_LOG_ASSERT(TrivialRowGrouping || dir.has_value(), "no optimization direction given!");
        if (!dir.has_value() || maximize(*dir)) {
            return RS<storm::OptimizationDirection::Maximize>(operand, offsets, numIterations, precision, iterationCallback);
        } else {
            return RS<storm::OptimizationDirection::Minimize>(operand, offsets, numIterations, precision, iterationCallback);
        }
    }

    auto RS(std::vector<TargetValueType>& operand, std::vector<TargetValueType> const& offsets, TargetValueType const& precision,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        uint64_t numIterations = 0;
        return RS(operand, offsets, numIterations, precision, dir, iterationCallback);
    }

   private:
    std::shared_ptr<ValueIterationOperator<ExactValueType, TrivialRowGrouping>> _exactOperator;
    std::shared_ptr<ValueIterationOperator<ImpreciseValueType, TrivialRowGrouping>> _impreciseOperator;
};

}  // namespace solver::helper
}  // namespace storm

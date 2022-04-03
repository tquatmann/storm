
#pragma once
#include <functional>
#include <optional>
#include <vector>

#include "storm/solver/SolverStatus.h"
#include "storm/solver/helper/ValueIterationOperator.h"
#include "storm/utility/Extremum.h"

namespace storm {
class Environment;

namespace solver::helper {

template<typename ValueType, bool Relative>
struct StandardVITerminationCriterion {
    ValueType const precision;
    bool operator()(ValueType const& oldValue, ValueType const& newValue) const {
        if constexpr (Relative) {
            return storm::utility::abs<ValueType>(oldValue - newValue) <= storm::utility::abs<ValueType>(precision * oldValue);
        } else {
            return storm::utility::abs<ValueType>(oldValue - newValue) <= precision;
        }
    }
};

template<typename ValueType, storm::OptimizationDirection Dir, class TermCrit>
class VIOperatorBackend {
   public:
    template<typename... TermCritArgs>
    VIOperatorBackend(TermCritArgs... termCritArgs) : _termCrit{termCritArgs...} {
        // intentionally empty
    }
    void startNewIteration() {
        _converged = true;
    }

    void firstRow(ValueType&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
        _best = std::move(value);
    }

    void nextRow(ValueType&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
        _best &= value;
    }

    void applyUpdate(ValueType& currValue, [[maybe_unused]] uint64_t rowGroup) {
        if (_converged) {
            _converged = _termCrit(currValue, *_best);
        }
        currValue = std::move(*_best);
    }

    void endOfIteration() const {
        // intentionally left empty.
    }

    bool converged() const {
        return _converged;
    }

    bool constexpr abort() const {
        return false;
    }

   private:
    storm::utility::Extremum<Dir, ValueType> _best;
    TermCrit const _termCrit;
    bool _converged{true};
};

template<typename ValueType, bool TrivialRowGrouping = false>
class ValueIterationHelper {
   public:
    ValueIterationHelper(std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> viOperator) : _operator(viOperator) {
        // Intentionally left empty
    }

    template<typename BackEndType>
    auto VI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, BackEndType&& backend,
            std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        SolverStatus status{SolverStatus::InProgress};
        while (status == SolverStatus::InProgress) {
            ++numIterations;
            if (_operator->template applyInPlace(operand, offsets, backend)) {
                status = SolverStatus::Converged;
            } else if (iterationCallback) {
                status = iterationCallback(status);
            }
        }
        return status;
    }

    template<storm::OptimizationDirection Dir, bool Relative>
    auto VI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, ValueType const& precision,
            std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        VIOperatorBackend<ValueType, Dir, StandardVITerminationCriterion<ValueType, Relative>> backend{precision};
        return VI(operand, offsets, numIterations, std::move(backend), iterationCallback);
    }

    template<storm::OptimizationDirection Dir>
    auto VI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative, ValueType const& precision,
            std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        if (relative) {
            return VI<Dir, true>(operand, offsets, numIterations, precision, iterationCallback);
        } else {
            return VI<Dir, false>(operand, offsets, numIterations, precision, iterationCallback);
        }
    }

    auto VI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative, ValueType const& precision,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        STORM_LOG_ASSERT(TrivialRowGrouping || dir.has_value(), "no optimization direction given!");
        if (!dir.has_value() || maximize(*dir)) {
            return VI<storm::OptimizationDirection::Maximize>(operand, offsets, numIterations, relative, precision, iterationCallback);
        } else {
            return VI<storm::OptimizationDirection::Minimize>(operand, offsets, numIterations, relative, precision, iterationCallback);
        }
    }

    auto VI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, bool relative, ValueType const& precision,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(SolverStatus const&)> const& iterationCallback = {}) {
        uint64_t numIterations = 0;
        return VI(operand, offsets, numIterations, relative, precision, dir, iterationCallback);
    }

   private:
    std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> _operator;
};

}  // namespace solver::helper
}  // namespace storm

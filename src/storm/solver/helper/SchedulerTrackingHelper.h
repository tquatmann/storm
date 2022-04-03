
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

template<typename ValueType, storm::OptimizationDirection Dir>
class SchedulerTrackingBackend {
   public:
    SchedulerTrackingBackend(std::vector<uint64_t>& schedulerStorage, std::vector<storm::storage::SparseMatrixIndexType> const& rowGroupIndices,
                             bool applyUpdates)
        : _schedulerStorage(schedulerStorage), _applyUpdates(applyUpdates), _rowGroupIndices(rowGroupIndices) {
        // intentionally empty
    }
    void startNewIteration() {
        _converged = true;
    }

    void firstRow(ValueType&& value, uint64_t rowGroup, uint64_t row) {
        _currChoice = row - _rowGroupIndices[rowGroup];
        _best = std::move(value);
    }

    void nextRow(ValueType&& value, uint64_t rowGroup, uint64_t row) {
        if (_best &= value) {
            _currChoice = row - _rowGroupIndices[rowGroup];
        } else if (*_best == value) {
            // For rows that are equally good, we prefer the previously selected one.
            // This is necessary, e.g., for policy iteration correctness.
            if (uint64_t rowChoice = row - _rowGroupIndices[rowGroup]; rowChoice == _schedulerStorage[rowGroup]) {
                _currChoice = rowChoice;
            }
        }
    }

    void applyUpdate(ValueType& currValue, uint64_t rowGroup) {
        if (_applyUpdates) {
            currValue = std::move(*_best);
        }
        auto& choice = _schedulerStorage[rowGroup];
        if (_converged) {
            _converged = choice == _currChoice;
        }
        choice = _currChoice;
    }

    void endOfIteration() const {}

    bool converged() const {
        return _converged;
    }

    bool constexpr abort() const {
        return false;
    }

   private:
    std::vector<uint64_t>& _schedulerStorage;
    bool const _applyUpdates;
    std::vector<storm::storage::SparseMatrixIndexType> const& _rowGroupIndices;

    bool _converged;
    storm::utility::Extremum<Dir, ValueType> _best;
    uint64_t _currChoice;
};

template<typename ValueType>
class SchedulerTrackingHelper {
   public:
    SchedulerTrackingHelper(std::shared_ptr<ValueIterationOperator<ValueType, false>> viOperator) : _operator(viOperator) {
        // Intentionally left empty
    }

    template<storm::OptimizationDirection Dir>
    bool computeScheduler(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, std::vector<uint64_t>& schedulerStorage, bool applyUpdates) {
        SchedulerTrackingBackend<ValueType, Dir> backend(schedulerStorage, _operator->getRowGroupIndices(), applyUpdates);
        return _operator->template applyInPlace(operand, offsets, backend);
    }

    bool computeScheduler(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, storm::OptimizationDirection const& dir,
                          std::vector<uint64_t>& schedulerStorage, bool applyUpdates) {
        if (maximize(dir)) {
            return computeScheduler<storm::OptimizationDirection::Maximize>(operand, offsets, schedulerStorage, applyUpdates);
        } else {
            return computeScheduler<storm::OptimizationDirection::Minimize>(operand, offsets, schedulerStorage, applyUpdates);
        }
    }

   private:
    std::shared_ptr<ValueIterationOperator<ValueType, false>> _operator;
};

}  // namespace solver::helper
}  // namespace storm

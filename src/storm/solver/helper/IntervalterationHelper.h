#pragma once
#include <type_traits>
#include <vector>
#include "storm/solver/SolverStatus.h"
#include "storm/solver/helper/ValueIterationHelper.h"
#include "storm/solver/helper/ValueIterationOperator.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/constants.h"
#include "storm/utility/vector.h"

namespace storm {
class Environment;

namespace solver::helper {

template<typename ValueType>
struct IIData {
    std::vector<ValueType> const &x, y;
    SolverStatus const status;
};

template<typename ValueType, bool TrivialRowGrouping = false>
class IntervalIterationHelper {
   public:
    IntervalIterationHelper(std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> viOperator) : _operator(viOperator) {
        // Intentionally left empty.
    }

    template<OptimizationDirection Dir>
    auto II(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative,
            ValueType const& precision, std::function<SolverStatus(IIData<ValueType> const&)> const& iterationCallback = {},
            std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        SolverStatus status{SolverStatus::InProgress};
        IIBackend<Dir> backend;
        uint64_t convergenceCheckState = 0;
        std::function<void()> getNextConvergenceCheckState;
        if (relevantValues) {
            convergenceCheckState = relevantValues->getNextSetIndex(0);
            getNextConvergenceCheckState = [&convergenceCheckState, &relevantValues]() {
                convergenceCheckState = relevantValues->getNextSetIndex(++convergenceCheckState);
            };
        } else {
            getNextConvergenceCheckState = [&convergenceCheckState]() { ++convergenceCheckState; };
        }
        while (status == SolverStatus::InProgress) {
            ++numIterations;
            _operator->template applyInPlace(xy, offsets, backend);
            if (checkConvergence(xy, convergenceCheckState, getNextConvergenceCheckState, relative, precision)) {
                status = SolverStatus::Converged;
            } else if (iterationCallback) {
                status = iterationCallback(IIData<ValueType>({xy.first, xy.second, status}));
            }
        }
        return status;
    }

    auto II(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative,
            ValueType const& precision, std::optional<storm::OptimizationDirection> const& dir,
            std::function<SolverStatus(IIData<ValueType> const&)> const& iterationCallback,
            std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        if (!dir.has_value() || maximize(*dir)) {
            return II<OptimizationDirection::Maximize>(xy, offsets, numIterations, relative, precision, iterationCallback, relevantValues);
        } else {
            return II<OptimizationDirection::Minimize>(xy, offsets, numIterations, relative, precision, iterationCallback, relevantValues);
        }
    }

    auto II(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative, ValueType const& precision,
            std::function<void(std::vector<ValueType>&)> const& prepareLowerBounds, std::function<void(std::vector<ValueType>&)> const& prepareUpperBounds,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(IIData<ValueType> const&)> const& iterationCallback = {},
            std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        // Create two vectors x and y using the given operand plus an auxiliary vector.
        std::pair<std::vector<ValueType>, std::vector<ValueType>> xy;
        auto& auxVector = _operator->allocateAuxiliaryVector(operand.size());
        xy.first.swap(operand);
        xy.second.swap(auxVector);
        prepareLowerBounds(xy.first);
        prepareUpperBounds(xy.second);
        auto doublePrec = precision + precision;
        if constexpr (std::is_same_v<ValueType, double>) {
            doublePrec -= precision * 1e-6;  // be slightly more precise to avoid a good chunk of floating point issues
        }
        auto status = II(xy, offsets, numIterations, relative, doublePrec, dir, iterationCallback, relevantValues);
        auto two = storm::utility::convertNumber<ValueType>(2.0);
        // get the average of lower- and upper result
        storm::utility::vector::applyPointwise<ValueType, ValueType, ValueType>(
            xy.first, xy.second, xy.first, [&two](ValueType const& a, ValueType const& b) -> ValueType { return (a + b) / two; });
        // Swap operand and aux vector back to original positions.
        xy.first.swap(operand);
        xy.second.swap(auxVector);
        _operator->freeAuxiliaryVector();
        return status;
    }

    auto II(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, bool relative, ValueType const& precision,
            std::function<void(std::vector<ValueType>&)> const& prepareLowerBounds, std::function<void(std::vector<ValueType>&)> const& prepareUpperBounds,
            std::optional<storm::OptimizationDirection> const& dir = {}, std::function<SolverStatus(IIData<ValueType> const&)> const& iterationCallback = {},
            std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        uint64_t numIterations = 0;
        return II(operand, offsets, numIterations, relative, precision, prepareLowerBounds, prepareUpperBounds, dir, iterationCallback, relevantValues);
    }

   private:
    template<OptimizationDirection Dir>
    class IIBackend {
       public:
        void startNewIteration() {}

        void firstRow(std::pair<ValueType, ValueType>&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            _xBest = std::move(value.first);
            _yBest = std::move(value.second);
        }

        void nextRow(std::pair<ValueType, ValueType>&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            assert(!TrivialRowGrouping);
            _xBest &= std::move(value.first);
            _yBest &= std::move(value.second);
        }

        void applyUpdate(ValueType& xCurr, ValueType& yCurr, [[maybe_unused]] uint64_t rowGroup) {
            xCurr = std::max(xCurr, *_xBest);
            yCurr = std::min(yCurr, *_yBest);
        }

        void endOfIteration() const {
            // intentionally left empty.
        }

        bool constexpr converged() const {
            return false;
        }

        bool constexpr abort() const {
            return false;
        }

       private:
        storm::utility::Extremum<Dir, ValueType> _xBest, _yBest;
    };

    bool checkConvergence(std::pair<std::vector<ValueType>, std::vector<ValueType>> const& xy, uint64_t& convergenceCheckState,
                          std::function<void()> const& getNextConvergenceCheckState, bool relative, ValueType const& precision) {
        if (relative) {
            for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                ValueType const& l = xy.first[convergenceCheckState];
                ValueType const& u = xy.second[convergenceCheckState];
                if (l > storm::utility::zero<ValueType>()) {
                    if ((u - l) > l * precision) {
                        return false;
                    }
                } else if (u < storm::utility::zero<ValueType>()) {
                    if ((l - u) < u * precision) {
                        return false;
                    }
                } else {  //  l <= 0 <= u
                    if (l != u) {
                        return false;
                    }
                }
            }
        } else {
            for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                if (xy.second[convergenceCheckState] - xy.first[convergenceCheckState] > precision) {
                    return false;
                }
            }
        }
        return true;
    }

    std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> _operator;
};

}  // namespace solver::helper
}  // namespace storm
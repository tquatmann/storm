#pragma once
#include <optional>
#include <type_traits>
#include <vector>

#include "storm/exceptions/UnexpectedException.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/ModelCheckerSettings.h"
#include "storm/solver/SolverStatus.h"
#include "storm/solver/helper/ValueIterationOperator.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/vector.h"

namespace storm {
class Environment;

namespace solver::helper {

template<typename ValueType, bool TrivialRowGrouping = false>
class SoundValueIterationHelper {
   public:
    SoundValueIterationHelper(std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> viOperator) : _operator(viOperator) {
        _sizeOfLargestRowGroup = 1;
        if constexpr (!TrivialRowGrouping) {
            auto it = _operator->getRowGroupIndices().cbegin();
            auto itEnd = _operator->getRowGroupIndices().cend() - 1;
            while (it != itEnd) {
                auto const curr = *it;
                _sizeOfLargestRowGroup = std::max(_sizeOfLargestRowGroup, *(++it) - curr);
            }
        }
    }

    enum class SVIStage { Initial, y_less_1, b_eq_d };

    template<OptimizationDirection Dir, SVIStage Stage>
    class SVIBackend {
       public:
        static const SVIStage CurrentStage = Stage;
        using RowValueStorageType = std::vector<std::pair<ValueType, ValueType>>;

        SVIBackend(RowValueStorageType& rowValueStorage, std::optional<ValueType> const& a, std::optional<ValueType> const& b,
                   std::optional<ValueType> const& d = {})
            : _currRowValues(rowValueStorage) {
            if (a.has_value()) {
                _a &= *a;
            }
            if (b.has_value()) {
                _b &= *b;
            }
            if (d.has_value()) {
                _d &= *d;
            }
        }

        void startNewIteration() {
            _allYLessOne = true;
            _curr_a.reset();
            _curr_b.reset();
        }

        void firstRow(std::pair<ValueType, ValueType>&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            assert(_currRowValuesIndex == 0);
            if constexpr (!TrivialRowGrouping) {
                _bestValue.reset();
            }
            _best = std::move(value);
        }

        void nextRow(std::pair<ValueType, ValueType>&& value, [[maybe_unused]] uint64_t rowGroup, [[maybe_unused]] uint64_t row) {
            assert(!TrivialRowGrouping);
            assert(_currRowValuesIndex < _currRowValues.size());
            if (Stage == SVIStage::Initial && _b.empty()) {
                if (value.second > _best.second || (value.second == _best.second && better(value.first, _best.first))) {
                    std::swap(value, _best);
                }
                _currRowValues[_currRowValuesIndex++] = std::move(value);
            } else {
                assert(!_b.empty());
                auto const& b = Stage == SVIStage::b_eq_d ? *_d : *_b;
                if (_bestValue.empty()) {
                    _bestValue = _best.first + b * _best.second;
                }
                if (ValueType currentValue = value.first + b * value.second; _bestValue &= currentValue) {
                    std::swap(value, _best);
                    if (Stage != SVIStage::b_eq_d && value.second < _best.second) {
                        // We need to store the 'old' best values as they might be relevant for the decision value.
                        _currRowValues[_currRowValuesIndex++] = std::move(value);
                    }
                } else if (_best.second > value.second) {
                    if (*_bestValue == currentValue) {
                        // In this case we have the same value, but the current row is to be preferred as it has a smaller y value
                        std::swap(value, _best);
                    } else if (Stage != SVIStage::b_eq_d) {
                        // In this case we have a worse weighted value
                        // However, this could be relevant for the decision value
                        _currRowValues[_currRowValuesIndex++] = std::move(value);
                    }
                }
            }
        }

        void applyUpdate(ValueType& xCurr, ValueType& yCurr, [[maybe_unused]] uint64_t rowGroup) {
            std::swap(xCurr, _best.first);
            std::swap(yCurr, _best.second);
            if constexpr (Stage != SVIStage::b_eq_d && !TrivialRowGrouping) {
                // Update decision value
                while (_currRowValuesIndex) {
                    if (auto const& rowVal = _currRowValues[--_currRowValuesIndex]; yCurr > rowVal.second) {
                        _d &= (rowVal.first - xCurr) / (yCurr - rowVal.second);
                    }
                }
            } else {
                assert(_currRowValuesIndex == 0);
            }

            // keep track of bounds a,b
            if constexpr (Stage == SVIStage::Initial) {
                if (_allYLessOne) {
                    if (yCurr < storm::utility::one<ValueType>()) {
                        ValueType val = xCurr / (storm::utility::one<ValueType>() - yCurr);
                        _curr_a &= val;
                        _curr_b &= val;
                    } else {
                        _allYLessOne = false;
                    }
                }
            } else {
                STORM_LOG_ASSERT(yCurr < storm::utility::one<ValueType>(), "Unexpected y value for this stage.");
                ValueType val = xCurr / (storm::utility::one<ValueType>() - yCurr);
                _curr_a &= val;
                _curr_b &= val;
            }
        }

        void endOfIteration() {
            _nextStage = Stage;
            if (_nextStage == SVIStage::Initial && _allYLessOne) {
                _nextStage = SVIStage::y_less_1;
            }
            if (_nextStage == SVIStage::y_less_1 || _nextStage == SVIStage::b_eq_d) {
                _a &= std::move(*_curr_a);
                if (_nextStage == SVIStage::y_less_1) {
                    _curr_b &= _d;
                    _b &= std::move(*_curr_b);
                    if (!_d.empty() && *_b == *_d) {
                        _nextStage = SVIStage::b_eq_d;
                    }
                } else {
                    // in the b_eq_d stage, we slightly repurpose _b and _d:
                    // _b is now used to track an upper bound (which can pass _d)
                    // _d is now used for the weighting when selecting the best row
                    _b &= std::move(*_curr_b);
                }
            }
        }

        bool constexpr converged() const {
            return false;
        }

        bool constexpr abort() const {
            return false;
        }

        std::optional<ValueType> a() const {
            return _a.getOptionalValue();
        }

        std::optional<ValueType> b() const {
            return _b.getOptionalValue();
        }

        std::optional<ValueType> d() const {
            return _d.getOptionalValue();
        }

        bool moveToNextStage() const {
            return _nextStage != Stage;
        }

        template<SVIStage NewStage>
        auto createBackendForNextStage() const {
            std::optional<ValueType> d;
            if (NewStage == SVIStage::b_eq_d && !_b.empty())
                d = *_b;
            else if (NewStage != SVIStage::Initial && !_d.empty())
                d = *_d;
            return SVIBackend<Dir, NewStage>(_currRowValues, a(), b(), d);
        }

        SVIStage const& nextStage() const {
            return _nextStage;
        }

       private:
        static bool better(ValueType const& lhs, ValueType const& rhs) {
            if constexpr (minimize(Dir)) {
                return lhs < rhs;
            } else {
                return lhs > rhs;
            }
        }

        using ExtremumDir = storm::utility::Extremum<Dir, ValueType>;
        using ExtremumInvDir = storm::utility::Extremum<invert(Dir), ValueType>;

        ExtremumDir _a, _d;
        ExtremumInvDir _b;

        SVIStage _nextStage{Stage};

        ExtremumDir _curr_b;
        ExtremumInvDir _curr_a;
        bool _allYLessOne;

        std::pair<ValueType, ValueType> _best;
        ExtremumDir _bestValue;
        RowValueStorageType& _currRowValues;
        uint64_t _currRowValuesIndex{0};
    };

    struct SVIData {
        SolverStatus status;
        std::pair<std::vector<ValueType>, std::vector<ValueType>> const& xy;
        std::optional<ValueType> const a, b;

        void trySetAverage(std::vector<ValueType>& out) const {
            if (a.has_value() && b.has_value()) {
                ValueType abAvg = (*a + *b) / storm::utility::convertNumber<ValueType, uint64_t>(2);
                storm::utility::vector::applyPointwise(xy.first, xy.second, out,
                                                       [&abAvg](ValueType const& xVal, ValueType const& yVal) -> ValueType { return xVal + abAvg * yVal; });
            }
        }

        void trySetLowerUpper(std::vector<ValueType>& lowerOut, std::vector<ValueType>& upperOut) const {
            auto [min, max] = std::minmax(*a, *b);
            uint64_t const size = xy.first.size();
            for (uint64_t i = 0; i < size; ++i) {
                // We allow setting both vectors "in-place", e.g. we might have &lowerOut == &xy.first.
                // This requires to use temporary values.
                ValueType xi = xy.first[i];
                ValueType yi = xy.second[i];
                lowerOut[i] = xi + min * yi;
                upperOut[i] = xi + max * yi;
            }
        }

        bool checkCustomTerminationCondition(storm::solver::TerminationCondition<ValueType> const& condition) const {
            if (a.has_value() && b.has_value()) {
                if (condition.requiresGuarantee(storm::solver::SolverGuarantee::GreaterOrEqual)) {
                    auto max = std::max(*a, *b);
                    return condition.terminateNow([&](uint64_t const& i) { return xy.first[i] + xy.second[i] * max; },
                                                  storm::solver::SolverGuarantee::GreaterOrEqual);
                } else if (condition.requiresGuarantee(storm::solver::SolverGuarantee::LessOrEqual)) {
                    auto min = std::min(*a, *b);
                    return condition.terminateNow([&](uint64_t const& i) { return xy.first[i] + xy.second[i] * min; },
                                                  storm::solver::SolverGuarantee::GreaterOrEqual);
                }
            }
            return false;
        }

        bool checkConvergence(uint64_t& convergenceCheckState, std::function<void()> const& getNextConvergenceCheckState, bool relative,
                              ValueType const& precision) {
            if (!a.has_value() || !b.has_value())
                return false;
            if (*a == *b)
                return true;
            if (relative) {
                auto [min, max] = std::minmax(*a, *b);
                if (min >= storm::utility::zero<ValueType>()) {
                    ValueType const val = (max - min) / precision - min;
                    for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                        if (!storm::utility::isZero(xy.second[convergenceCheckState]) &&
                            val > xy.first[convergenceCheckState] / xy.second[convergenceCheckState]) {
                            return false;
                        }
                    }
                } else if (max <= storm::utility::zero<ValueType>()) {
                    ValueType const val = (min - max) / precision - max;
                    for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                        if (!storm::utility::isZero(xy.second[convergenceCheckState]) &&
                            val < xy.first[convergenceCheckState] / xy.second[convergenceCheckState]) {
                            return false;
                        }
                    }
                } else {
                    for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                        ValueType l = xy.first[convergenceCheckState] + min * xy.second[convergenceCheckState];
                        ValueType u = xy.first[convergenceCheckState] + max * xy.second[convergenceCheckState];
                        assert(u >= l);
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
                }
            } else {
                ValueType val = precision / storm::utility::abs<ValueType>(*b - *a);
                for (; convergenceCheckState < xy.first.size(); getNextConvergenceCheckState()) {
                    if (xy.second[convergenceCheckState] > val) {
                        return false;
                    }
                }
            }
            return true;
        }
    };

    template<typename BackendType>
    auto SVI(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::pair<std::vector<ValueType> const*, ValueType> const& offsets,
             uint64_t& numIterations, bool relative, ValueType const& precision, BackendType&& backend,
             std::function<SolverStatus(SVIData const&)> const& iterationCallback, std::optional<storm::storage::BitVector> const& relevantValues,
             uint64_t convergenceCheckState = 0) {
        if constexpr (BackendType::CurrentStage == SVIStage::Initial) {
            xy.first.assign(xy.first.size(), storm::utility::zero<ValueType>());
            xy.second.assign(xy.first.size(), storm::utility::one<ValueType>());
            convergenceCheckState = relevantValues.has_value() ? relevantValues->getNextSetIndex(0ull) : 0ull;
        }
        std::function<void()> getNextConvergenceCheckState;
        if (relevantValues) {
            getNextConvergenceCheckState = [&convergenceCheckState, &relevantValues]() {
                convergenceCheckState = relevantValues->getNextSetIndex(++convergenceCheckState);
            };
        } else {
            getNextConvergenceCheckState = [&convergenceCheckState]() { ++convergenceCheckState; };
        }

        while (true) {
            ++numIterations;
            _operator->template applyInPlace(xy, offsets, backend);
            SVIData data{SolverStatus::InProgress, xy, backend.a(), backend.b()};
            if (data.checkConvergence(convergenceCheckState, getNextConvergenceCheckState, relative, precision)) {
                return SVIData{SolverStatus::Converged, xy, backend.a(), backend.b()};
            } else {
                if (iterationCallback) {
                    SVIData data{SolverStatus::InProgress, xy, backend.a(), backend.b()};
                    data.status = iterationCallback(data);
                    if (data.status != SolverStatus::InProgress) {
                        return data;
                    }
                }
                if (backend.moveToNextStage()) {
                    switch (backend.nextStage()) {
                        case SVIStage::y_less_1:
                            return SVI(xy, offsets, numIterations, relative, precision, backend.template createBackendForNextStage<SVIStage::y_less_1>(),
                                       iterationCallback, relevantValues, convergenceCheckState);
                        case SVIStage::b_eq_d:
                            return SVI(xy, offsets, numIterations, relative, precision, backend.template createBackendForNextStage<SVIStage::b_eq_d>(),
                                       iterationCallback, relevantValues, convergenceCheckState);
                        default:
                            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected next stage");
                    }
                }
            }
        }
    }

    template<storm::OptimizationDirection Dir>
    auto SVI(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::pair<std::vector<ValueType> const*, ValueType> const& offsets,
             uint64_t& numIterations, bool relative, ValueType const& precision, std::optional<ValueType> const& a, std::optional<ValueType> const& b,
             std::function<SolverStatus(SVIData const&)> const& iterationCallback = {}, std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        typename SVIBackend<Dir, SVIStage::Initial>::RowValueStorageType rowValueStorage;
        rowValueStorage.resize(_sizeOfLargestRowGroup - 1);
        return SVI(xy, offsets, numIterations, relative, precision, SVIBackend<Dir, SVIStage::Initial>(rowValueStorage, a, b), iterationCallback,
                   relevantValues);
    }

    auto SVI(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::pair<std::vector<ValueType> const*, ValueType> const& offsets,
             uint64_t& numIterations, bool relative, ValueType const& precision, std::optional<storm::OptimizationDirection> const& dir,
             std::optional<ValueType> const& lowerBound, std::optional<ValueType> const& upperBound,
             std::function<SolverStatus(SVIData const&)> const& iterationCallback, std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        if (!dir.has_value() || maximize(*dir)) {
            // When we maximize, a is the lower bound and b is the upper bound
            return SVI<storm::OptimizationDirection::Maximize>(xy, offsets, numIterations, relative, precision, lowerBound, upperBound, iterationCallback,
                                                               relevantValues);
        } else {
            // When we minimize, b is the lower bound and a is the upper bound
            return SVI<storm::OptimizationDirection::Minimize>(xy, offsets, numIterations, relative, precision, upperBound, lowerBound, iterationCallback,
                                                               relevantValues);
        }
    }

    auto SVI(std::pair<std::vector<ValueType>, std::vector<ValueType>>& xy, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative,
             ValueType const& precision, std::optional<storm::OptimizationDirection> const& dir, std::optional<ValueType> const& lowerBound,
             std::optional<ValueType> const& upperBound, std::function<SolverStatus(SVIData const&)> const& iterationCallback,
             std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        std::pair<std::vector<ValueType> const*, ValueType> offsetsPair{&offsets, storm::utility::zero<ValueType>()};
        return SVI(xy, offsetsPair, numIterations, relative, precision, dir, lowerBound, upperBound, iterationCallback, relevantValues);
    }

    auto SVI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, uint64_t& numIterations, bool relative, ValueType const& precision,
             std::optional<storm::OptimizationDirection> const& dir = {}, std::optional<ValueType> const& lowerBound = {},
             std::optional<ValueType> const& upperBound = {}, std::function<SolverStatus(SVIData const&)> const& iterationCallback = {},
             std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        // Create two vectors x and y using the given operand plus an auxiliary vector.
        std::pair<std::vector<ValueType>, std::vector<ValueType>> xy;
        auto& auxVector = _operator->allocateAuxiliaryVector(operand.size());
        xy.first.swap(operand);
        xy.second.swap(auxVector);
        auto doublePrec = precision + precision;
        if constexpr (std::is_same_v<ValueType, double>) {
            doublePrec -= precision * 1e-6;  // be slightly more precise to avoid a good chunk of floating point issues
        }
        auto res = SVI(xy, offsets, numIterations, relative, doublePrec, dir, lowerBound, upperBound, iterationCallback, relevantValues);
        res.trySetAverage(xy.first);
        // Swap operand and aux vector back to original positions.
        xy.first.swap(operand);
        xy.second.swap(auxVector);
        _operator->freeAuxiliaryVector();
        return res.status;
    }

    auto SVI(std::vector<ValueType>& operand, std::vector<ValueType> const& offsets, bool relative, ValueType const& precision,
             std::optional<storm::OptimizationDirection> const& dir = {}, std::optional<ValueType> const& lowerBound = {},
             std::optional<ValueType> const& upperBound = {}, std::function<SolverStatus(SVIData const&)> const& iterationCallback = {},
             std::optional<storm::storage::BitVector> const& relevantValues = {}) {
        uint64_t numIterations = 0;
        return SVI(operand, offsets, numIterations, relative, precision, dir, lowerBound, upperBound, iterationCallback, relevantValues);
    }

    std::shared_ptr<ValueIterationOperator<ValueType, TrivialRowGrouping>> _operator;
    uint64_t _sizeOfLargestRowGroup;
};

}  // namespace solver::helper
}  // namespace storm
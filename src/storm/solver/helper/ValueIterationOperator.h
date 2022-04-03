#pragma once
#include <vector>

#include "storm/storage/SparseMatrix.h"
#include "storm/utility/macros.h"

namespace storm {
class Environment;

namespace solver::helper {

template<typename ValueType, bool TrivialRowGrouping = false>
class ValueIterationOperator {
   public:
    using IndexType = uint64_t;

    void setMatrixBackwards(storm::storage::SparseMatrix<ValueType> const& matrix,
                            std::vector<storm::storage::SparseMatrixIndexType> const* rowGroupIndices = nullptr) {
        if constexpr (TrivialRowGrouping) {
            STORM_LOG_ASSERT(matrix.hasTrivialRowGrouping(), "Expected a matrix with trivial row grouping");
            STORM_LOG_ASSERT(rowGroupIndices == nullptr, "Row groups given, but grouping is supposed to be trivial.");
            _rowGroupIndices = nullptr;
        } else {
            if (rowGroupIndices) {
                _rowGroupIndices = rowGroupIndices;
            } else {
                _rowGroupIndices = &matrix.getRowGroupIndices();
            }
        }
        auto numRows = matrix.getRowCount();
        _matrixValues.clear();
        _matrixColumns.clear();
        _matrixValues.reserve(matrix.getNonzeroEntryCount());
        _matrixColumns.reserve(matrix.getNonzeroEntryCount() + numRows + 1);  // _matrixColumns also contain indications for when a row(group) starts
        if constexpr (!TrivialRowGrouping) {
            auto groupIt = _rowGroupIndices->crbegin();
            _matrixColumns.push_back(StartOfRowGroupIndicator);  // indicate start of first row(group)
            for (IndexType rowGroupEnd = *groupIt; rowGroupEnd != 0; rowGroupEnd = *groupIt) {
                IndexType rowIndex = *(++groupIt);
                STORM_LOG_ASSERT(rowIndex != rowGroupEnd, "There is an empty row group. This is not expected.");
                for (; rowIndex != rowGroupEnd; ++rowIndex) {
                    for (auto const& entry : matrix.getRow(rowIndex)) {
                        _matrixValues.push_back(entry.getValue());
                        _matrixColumns.push_back(entry.getColumn());
                    }
                    _matrixColumns.push_back(StartOfRowIndicator);  // Indicate start of next row
                }
                _matrixColumns.back() = StartOfRowGroupIndicator;  // Actually, this is the start of the next row group
            }
        } else {
            _matrixColumns.push_back(StartOfRowIndicator);  // Indicate start of first row
            for (IndexType rowIndex = numRows; rowIndex != 0;) {
                for (auto const& entry : matrix.getRow(--rowIndex)) {
                    _matrixValues.push_back(entry.getValue());
                    _matrixColumns.push_back(entry.getColumn());
                }
                _matrixColumns.push_back(StartOfRowIndicator);  // Indicate start of next row
            }
        }
    }

    template<typename OperandType, typename OffsetType, typename BackendType>
    bool applyInPlace(OperandType& operand, OffsetType const& offsets, BackendType& backend) {
        if (_hasSkippedRows) {
            return applyInPlace<OperandType, OffsetType, BackendType, true>(operand, offsets, backend);
        } else {
            return applyInPlace<OperandType, OffsetType, BackendType, false>(operand, offsets, backend);
        }
    }

    template<typename OperandType, typename OffsetType, typename BackendType, bool SkipIgnoredRows>
    bool applyInPlace(OperandType& operand, OffsetType const& offsets, BackendType& backend) {
        auto operandSize = getSize(operand);
        STORM_LOG_ASSERT(TrivialRowGrouping || _rowGroupIndices->size() == operandSize + 1, "Dimension mismatch");
        backend.startNewIteration();
        auto matrixValueIt = _matrixValues.cbegin();
        auto matrixColumnIt = _matrixColumns.cbegin();
        [[maybe_unused]] auto groupIt =
            TrivialRowGrouping ? typename std::vector<storm::storage::SparseMatrixIndexType>::const_reverse_iterator() : _rowGroupIndices->crbegin();
        for (uint64_t opIndex = operandSize; opIndex != 0;) {
            --opIndex;
            if constexpr (TrivialRowGrouping) {
                backend.firstRow(applyRow(matrixColumnIt, matrixValueIt, operand, offsets, opIndex), opIndex, opIndex);
            } else {
                auto offsetIndex = *(++groupIt);
                if constexpr (SkipIgnoredRows) {
                    offsetIndex += skipMultipleIgnoredRows(matrixColumnIt, matrixValueIt);
                }
                backend.firstRow(applyRow(matrixColumnIt, matrixValueIt, operand, offsets, offsetIndex), opIndex, offsetIndex);
                while (*matrixColumnIt < StartOfRowGroupIndicator) {
                    ++offsetIndex;
                    if (!SkipIgnoredRows || !skipIgnoredRow(matrixColumnIt, matrixValueIt)) {
                        backend.nextRow(applyRow(matrixColumnIt, matrixValueIt, operand, offsets, offsetIndex), opIndex, offsetIndex);
                    }
                }
            }
            if constexpr (isPair<OperandType>::value) {
                backend.applyUpdate(operand.first[opIndex], operand.second[opIndex], opIndex);
            } else {
                backend.applyUpdate(operand[opIndex], opIndex);
            }
            if (backend.abort()) {
                return backend.converged();
            }
        }
        STORM_LOG_ASSERT(matrixColumnIt + 1 == _matrixColumns.cend(), "Unexpected position of matrix column iterator.");
        STORM_LOG_ASSERT(matrixValueIt == _matrixValues.cend(), "Unexpected position of matrix column iterator.");
        backend.endOfIteration();
        return backend.converged();
    }

    void unsetIgnoredRows() {
        for (auto& c : _matrixColumns) {
            if (c >= StartOfRowIndicator) {
                c &= StartOfRowGroupIndicator;
            }
        }
        _hasSkippedRows = false;
    }

    template<typename IgnoredCallback>
    void setIgnoredRows(bool useLocalRowIndices, IgnoredCallback const& ignore) {
        static_assert(!TrivialRowGrouping);
        IndexType groupIndex{_rowGroupIndices->size() - 1};
        IndexType rowIndex{};
        auto colIt = _matrixColumns.begin();
        while (true) {
            STORM_LOG_ASSERT(colIt != _matrixColumns.end(), "VI Operator in invalid state.");
            STORM_LOG_ASSERT(*colIt >= StartOfRowIndicator, "VI Operator in invalid state.");
            if (*colIt >= StartOfRowGroupIndicator) {
                if (groupIndex == 0) {
                    break;
                }
                --groupIndex;
                rowIndex = useLocalRowIndices ? 0 : (*_rowGroupIndices)[groupIndex];
            } else {
                ++rowIndex;
            }

            if (!ignore(groupIndex, rowIndex)) {
                *colIt &= StartOfRowGroupIndicator;  // Clear number of skipped entries
                moveToEndOfRow(colIt);
            } else if ((*colIt & SkipNumEntriesMask) == 0) {  // i.e. should ignore but is not already ignored
                auto currColIt = colIt;
                moveToEndOfRow(colIt);
                *currColIt += std::distance(currColIt, colIt) - 1;  // set number of skipped entries
            }
        }
        _hasSkippedRows = false;
    }

    auto const& getRowGroupIndices() const {
        static_assert(!TrivialRowGrouping);
        return *_rowGroupIndices;
    }

    std::vector<ValueType>& allocateAuxiliaryVector(uint64_t size, std::optional<ValueType> const& initialValue = {}) {
        STORM_LOG_ASSERT(!_auxiliaryVectorUsedExternally, "Auxiliary vector already in use.");
        if (initialValue) {
            _auxiliaryVector.assign(size, *initialValue);
        } else {
            _auxiliaryVector.resize(size);
        }
        _auxiliaryVectorUsedExternally = true;
        return _auxiliaryVector;
    }

    void freeAuxiliaryVector() {
        _auxiliaryVectorUsedExternally = false;
    }

   private:
    template<typename OpT, typename OffT>
    OpT initializeRowRes(std::vector<OpT> const&, std::vector<OffT> const& offsets, uint64_t offsetIndex) {
        return offsets[offsetIndex];
    }

    template<typename OpT1, typename OpT2, typename OffT>
    std::pair<OpT1, OpT2> initializeRowRes(std::pair<std::vector<OpT1>, std::vector<OpT2>> const&, std::vector<OffT> const& offsets, uint64_t offsetIndex) {
        return {offsets[offsetIndex], offsets[offsetIndex]};
    }

    template<typename OpT1, typename OpT2, typename OffT1, typename OffT2>
    std::pair<OpT1, OpT2> initializeRowRes(std::pair<std::vector<OpT1>, std::vector<OpT2>> const&, std::pair<std::vector<OffT1> const*, OffT2> const& offsets,
                                           uint64_t offsetIndex) {
        return {(*offsets.first)[offsetIndex], offsets.second};
    }

    template<typename OperandType, typename OffsetType>
    auto applyRow(std::vector<IndexType>::const_iterator& matrixColumnIt, typename std::vector<ValueType>::const_iterator& matrixValueIt,
                  OperandType const& operand, OffsetType const& offsets, uint64_t offsetIndex) {
        STORM_LOG_ASSERT(*matrixColumnIt >= StartOfRowIndicator, "VI Operator in invalid state.");
        auto result{initializeRowRes(operand, offsets, offsetIndex)};
        for (++matrixColumnIt; *matrixColumnIt < StartOfRowIndicator; ++matrixColumnIt, ++matrixValueIt) {
            if constexpr (isPair<OperandType>::value) {
                result.first += operand.first[*matrixColumnIt] * (*matrixValueIt);
                result.second += operand.second[*matrixColumnIt] * (*matrixValueIt);
            } else {
                result += operand[*matrixColumnIt] * (*matrixValueIt);
            }
        }
        return result;
    }

    void moveToEndOfRow(std::vector<IndexType>::iterator& matrixColumnIt) {
        do {
            ++matrixColumnIt;
        } while (*matrixColumnIt < StartOfRowIndicator);
    }

    bool skipIgnoredRow(std::vector<IndexType>::const_iterator& matrixColumnIt, typename std::vector<ValueType>::const_iterator& matrixValueIt) {
        if (IndexType entriesToSkip = (*matrixColumnIt & SkipNumEntriesMask)) {
            matrixColumnIt += entriesToSkip + 1;
            matrixValueIt += entriesToSkip;
            return true;
        }
        return false;
    }

    uint64_t skipMultipleIgnoredRows(std::vector<IndexType>::const_iterator& matrixColumnIt, typename std::vector<ValueType>::const_iterator& matrixValueIt) {
        IndexType result{0ull};
        while (skipIgnoredRow(matrixColumnIt, matrixValueIt)) {
            ++result;
            STORM_LOG_ASSERT(*matrixColumnIt >= StartOfRowIndicator, "Undexpected state of VI operator");
            // We (currently) don't use this past the end of a row group, so we may have this additional sanity check:
            STORM_LOG_ASSERT(*matrixColumnIt < StartOfRowGroupIndicator, "Undexpected state of VI operator");
        }
        return result;
    }

    // Auxiliary helpers used for metaprogramming

    template<typename T>
    uint64_t getSize(std::vector<T> const& vec) {
        return vec.size();
    }

    template<typename T1, typename T2>
    uint64_t getSize(std::pair<T1, T2> const& pairOfVec) {
        return pairOfVec.first.size();
    }

    template<typename>
    struct isPair : std::false_type {};

    template<typename T1, typename T2>
    struct isPair<std::pair<T1, T2>> : std::true_type {};

    IndexType const StartOfRowIndicator = 1ull << 63;                               // 10000..0
    IndexType const StartOfRowGroupIndicator = StartOfRowIndicator + (1ull << 62);  // 11000..0
    IndexType const SkipNumEntriesMask = ~StartOfRowGroupIndicator;                 // 00111..1

    std::vector<ValueType> _matrixValues;
    std::vector<IndexType> _matrixColumns;
    std::vector<storm::storage::SparseMatrixIndexType> const* _rowGroupIndices;

    bool _hasSkippedRows{false};

    std::vector<ValueType> _auxiliaryVector;
    bool _auxiliaryVectorUsedExternally{false};
};

}  // namespace solver::helper
}  // namespace storm

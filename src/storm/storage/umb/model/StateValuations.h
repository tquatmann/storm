#pragma once

#include <bit>
#include <bitset>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <ranges>
#include <span>

#include "storm/storage/expressions/ExpressionManager.h"
#include "storm/storage/expressions/Variable.h"
#include "storm/storage/umb/model/ValuationDescription.h"
#include "storm/storage/umb/model/ValueEncoding.h"
#include "storm/utility/bitoperations.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"

namespace storm::umb {

class StateValuations {
   public:
    class Valuation {
        friend class StateValuations;  // Only StateValuations can construct Valuation

       public:
        void forEachAssignment(auto const& callback) const {
            DescriptionIterator descriptionIterator = description.variables.cbegin();
            BitIndex currentValueStart = 0;
            advancePadding(currentValueStart, descriptionIterator);
            while (descriptionIterator != description.variables.cend()) {
                invokeForCurrent<false>(currentValueStart, descriptionIterator, callback);
                ++descriptionIterator;
                advancePadding(currentValueStart, descriptionIterator);
                STORM_LOG_ASSERT(currentValueStart <= bytes.size() * 8, "Unexpected end of valuation data.");
            }
        }

       private:
        Valuation(std::span<char const> bytes, ValuationDescription const& description) : bytes(bytes), description(description) {}

        using BitIndex = uint64_t;
        using DescriptionIterator = typename decltype(ValuationDescription::variables)::const_iterator;
        using Integer = storm::NumberTraits<storm::RationalNumber>::IntegerType;

        void advancePadding(BitIndex& currentValueStart, DescriptionIterator& descriptionIterator) const {
            if (descriptionIterator != description.variables.cend() && std::holds_alternative<ValuationDescription::Padding>(*descriptionIterator)) {
                currentValueStart += std::get<ValuationDescription::Padding>(*descriptionIterator).padding;
                ++descriptionIterator;
                advancePadding(currentValueStart, descriptionIterator);
            }
        }

        uint64_t readUint64(BitIndex start, uint64_t const size) const {
            STORM_LOG_ASSERT(size <= 64, "Invalid bit range.");
            auto const firstByte = start / 8;
            auto const numBytes = (start % 8 + size + 7) / 8;
            uint64_t result;
            if (numBytes <= 8) {
                std::memcpy(&result, &bytes[firstByte], numBytes);
                result >>= (start % 8);
            } else {
                std::memcpy(&result, &bytes[firstByte], 8);
                result >>= (start % 8);
                uint64_t upperByte = std::bit_cast<uint8_t>(bytes[firstByte + 8]);
                result |= (upperByte << (64 - (start % 8)));
            }
            if (size < 64) {
                uint64_t const relevantBitMask = (1ull << size) - 1;
                result &= relevantBitMask;
            }
            return result;
        }

        template<bool Signed>
        Integer readInteger(BitIndex start, uint64_t const size) const {
            auto const num64BitChunks = (size + 63) / 64;
            auto chunksView = std::ranges::iota_view(0ull, num64BitChunks) | std::ranges::views::transform([&start, &size, this](auto i) -> uint64_t {
                                  return readUint64(start + i * 64, std::min<uint64_t>(64, size - i * 64));
                              });
            Integer result = ValueEncoding::decodeArbitraryPrecisionInteger<false>(chunksView);
            if constexpr (Signed) {
                // Check if this number is supposed to be negative
                if (result >= storm::utility::pow<Integer>(2, size - 1)) {
                    return result - storm::utility::pow<Integer>(2, size);
                }
            }
            return result;
        }

        bool safeAdd(int64_t a, int64_t b, int64_t& result) const {
            if (b > 0 && a > std::numeric_limits<int64_t>::max() - b) {
                return false;
            }
            if (b < 0 && a < std::numeric_limits<int64_t>::min() - b) {
                return false;
            }
            result = a + b;
            return true;
        }

        template<bool onlyTrivialTypes>
        void invokeForCurrent(BitIndex& currentValueStart, DescriptionIterator const& descriptionIterator, auto const& callback) const {
            using enum storm::umb::Type;
            STORM_LOG_ASSERT(std::holds_alternative<ValuationDescription::Variable>(*descriptionIterator), "Unexpected type of variable description.");
            auto const& variableDescription = std::get<ValuationDescription::Variable>(*descriptionIterator);
            STORM_LOG_ASSERT(currentValueStart < bytes.size() * 8, "Unexpected bit index.");
            auto const varBitSize = variableDescription.type.bitSize();
            if (variableDescription.type.bitSize() <= 64) {
                // Fast path for fixed-size, up to 64 bit values
                uint64_t rawContent = readUint64(currentValueStart, varBitSize);
                currentValueStart += varBitSize;
                switch (variableDescription.type.type) {
                    case Bool:
                        callback(rawContent != 0);
                        return;
                    case Uint:
                        if (std::in_range<int64_t>(rawContent)) {
                            int64_t result = rawContent;
                            if (safeAdd(result, variableDescription.offset.value_or(0), result)) {
                                callback(result);
                                return;
                            }
                        }
                        break;
                    case Int: {
                        uint64_t const mostSignificantBitMask = 1ull << (varBitSize - 1);
                        if (rawContent & mostSignificantBitMask) {
                            // Negative value: Two's complement (e.g. 1111...1101 is -3)
                            callback(-static_cast<int64_t>(~rawContent & (mostSignificantBitMask - 1)) - 1);
                            return;
                        } else {
                            // Positive value
                            callback(static_cast<int64_t>(rawContent));
                            return;
                        }
                    } break;
                    case Double:
                        callback(static_cast<double>(std::bit_cast<double>(rawContent) + variableDescription.offset.value_or(0)));
                        return;
                }
                // reaching this point means that we could not handle the value in the fast path
                // reset the current bit index to re-process the value in the slow path
                currentValueStart -= varBitSize;
            }
            // Handle larger or variable-size values
            if constexpr (onlyTrivialTypes) {
                STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                                "Non-trivial type or value of variable " << variableDescription.name << " is not supported in this context.");
            } else {
                switch (variableDescription.type.type) {
                    case Bool:
                        callback(readInteger<false>(currentValueStart, varBitSize) != Integer(0));
                        return;
                    case Uint:
                        callback(Integer(readInteger<false>(currentValueStart, varBitSize) +
                                         storm::utility::convertNumber<Integer>(variableDescription.offset.value_or(0))));
                        return;
                    case Int:
                        callback(Integer(readInteger<true>(currentValueStart, varBitSize) +
                                         storm::utility::convertNumber<Integer>(variableDescription.offset.value_or(0))));
                        return;
                    case Rational: {
                        auto const x = readUint64(currentValueStart, 16);
                        currentValueStart += 16;
                        storm::RationalNumber numerator = readInteger<true>(currentValueStart, x);
                        currentValueStart += x;
                        storm::RationalNumber denominator = readInteger<false>(currentValueStart, x);
                        currentValueStart += x;
                        callback(storm::RationalNumber(numerator / denominator));
                        return;
                    }
                    case String: {
                        STORM_LOG_ASSERT(currentValueStart % 8 == 0, "String values must be byte-aligned.");
                        auto const b = readUint64(currentValueStart, 16);
                        currentValueStart += 16;
                        auto substring = std::string_view(bytes.data() + currentValueStart / 8, b);
                        currentValueStart += b * 8;  // b is in bytes
                        callback(substring);
                        return;
                    }
                }
            }
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Type for variable " << variableDescription.name << " is not handled.");
        }

        std::span<char const> bytes;
        ValuationDescription const& description;
    };

    StateValuations(std::vector<ValuationDescription> description, std::vector<char> valuations, std::optional<std::vector<uint32_t>> classes = {},
                    std::optional<std::vector<uint64_t>> stringMapping = {}, std::optional<std::vector<char>> strings = {})
        : bytes(std::move(valuations)),
          description(std::move(description.front())),
          bitsize(this->description.sizeInBits()),
          expressionManager(std::make_shared<storm::expressions::ExpressionManager>()) {
        STORM_LOG_ASSERT(bitsize > 0, "Invalid alignment 0 for state valuations.");
        STORM_LOG_ASSERT(this->bytes.size() % bitsize == 0, "Corrupted state valuation data.");
        STORM_LOG_ASSERT(!this->stateValuationIndices || this->stateValuationIndices->back() == this->bytes.size() / bitsize,
                         "Corrupted state valuation index data: " << this->stateValuationIndices->back() << " vs. " << this->bytes.size() / bitsize);
        for (auto& varDesc : this->description.variables) {
            if (std::holds_alternative<ValuationDescription::Variable>(varDesc)) {
                variables.push_back(createVariable(*expressionManager, std::get<ValuationDescription::Variable>(varDesc)));
            }
        }
    }
    
    uint64_t getNumberOfStates() const {
        if (stateValuationIndices) {
            return stateValuationIndices->size() - 1;
        } else {
            return bytes.size() / bitsize;
        }
    }

    Valuation getValuation(uint64_t stateIndex) const {
        if (stateValuationIndices) {
            STORM_LOG_ASSERT(stateIndex < getNumberOfStates(),
                             "Invalid state index " << stateIndex << ". There are only " << getNumberOfStates() << " indexed states.");
            uint64_t const start = (*stateValuationIndices)[stateIndex] * bitsize;
            uint64_t const end = (*stateValuationIndices)[stateIndex + 1] * bitsize;
            STORM_LOG_ASSERT(end <= bytes.size(), "Corrupted state valuation data.");
            return Valuation(std::span<char const>(bytes.data() + start, end - start), description);
        } else {
            uint64_t const start = stateIndex * bitsize;
            STORM_LOG_ASSERT(bytes.size() <= start + bitsize, "Corrupted state valuation data.");
            return Valuation(std::span<char const>(bytes.data() + start, bitsize), description);
        }
    }

    auto getRange() {
        return std::ranges::iota_view(0ull, getNumberOfStates()) | std::ranges::views::transform([this](auto i) -> Valuation { return getValuation(i); });
    }

   private:
    struct Variable {
        storm::expressions::Variable const expressionVariable;
        std::reference_wrapper<ValuationDescription::Variable const> description;
    };

    static Variable createVariable(storm::expressions::ExpressionManager& expressionManager, ValuationDescription::Variable const& variableDescription) {
        storm::expressions::Type variableType;
        using enum storm::umb::Type;
        switch (variableDescription.type.type) {
            case Bool:
                variableType = expressionManager.getBooleanType();
                break;
            case Uint:
            case Int:
                variableType = expressionManager.getIntegerType();
                break;
            case Double:
            case Rational:
                variableType = expressionManager.getRationalType();
                break;
            case String:
                variableType = expressionManager.getStringType();
                break;
            default:
                STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected variable type.");
        }
        return {expressionManager.declareVariable(variableDescription.name, variableType), variableDescription};
    }

    std::vector<char> const bytes;
    std::optional<std::vector<uint64_t>> const stateValuationIndices;
    ValuationDescription const description;
    uint64_t const bitsize;
    std::shared_ptr<storm::expressions::ExpressionManager> expressionManager;

    std::vector<Variable> variables;
};
}  // namespace storm::umb
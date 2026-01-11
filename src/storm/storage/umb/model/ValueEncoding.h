#pragma once
#include <ranges>
#include <span>
#include <vector>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/storage/umb/model/GenericVector.h"
#include "storm/storage/umb/model/Types.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"

namespace storm::umb {

class ValueEncoding {
   public:
    /*!
     * returns func(<decoded_input>) where <decoded_input> is a range that (if necessary) decodes and converts the given input.
     * @param func a function that takes a range of value type and returns the result
     * @param input the input vector
     * @param sourceType The type that the input vector represents
     * @param csr a CSR that is used to decode the input vector. This is only necessary for non-trivially encoded source types (like rational).
     * @return
     */
    template<typename ValueType, typename SourceTypeEnum>
        requires std::is_enum_v<typename SourceTypeEnum::E>
    static auto applyDecodedVector(auto&& func, storm::umb::GenericVector const& input, SourceTypeEnum sourceType, CSR const& csr = {}) {
        using E = typename decltype(sourceType)::E;
        static_assert(std::is_enum_v<E>, "Source type must be an enum.");

        // find out how to interpret the input values
        switch (sourceType) {
            case E::Double:
                STORM_LOG_WARN_COND(
                    !storm::NumberTraits<ValueType>::IsExact,
                    "Some values are given in type double but will be converted to an exact (arbitrary precision) type. Rounding errors may occur.");
                STORM_LOG_ASSERT(input.template isType<double>(), "Unexpected type for values. Expected double.");
                return func(conversionView<ValueType>(input.template get<double>()));
            case E::Rational:
                STORM_LOG_WARN_COND(storm::NumberTraits<ValueType>::IsExact,
                                    "Some values are given in an exact type but converted to an inexact type. Rounding errors may occur.");
                if (input.template isType<storm::RationalNumber>()) {
                    return func(conversionView<ValueType>(input.template get<storm::RationalNumber>()));
                } else {
                    STORM_LOG_ASSERT(input.template isType<uint64_t>(), "Unexpected type for rational representation. Expected uint64.");
                    // Only this case might require the csr. It is ignored in all other cases.
                    if (csr) {
                        return func(conversionView<ValueType>(uint64ToRationalRangeView(input.template get<uint64_t>(), *csr)));
                    } else {
                        return func(conversionView<ValueType>(uint64ToRationalRangeView(input.template get<uint64_t>())));
                    }
                }
            case E::DoubleInterval:
                // For intervals, there is no suitable value conversion since we would drop the uncertainty
                if constexpr (!std::is_same_v<ValueType, storm::Interval>) {
                    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                                    "Some values are given as double intervals but a model with a non-interval type is requested.");
                    return func(std::ranges::empty_view<ValueType>{});
                } else {
                    if (input.template isType<storm::Interval>()) {
                        return func(input.template get<storm::Interval>());
                    } else {
                        STORM_LOG_ASSERT(input.template isType<double>(), "Unexpected type for double interval representation. Expected double.");
                        return func(doubleToIntervalRangeView(input.template get<double>()));
                    }
                }
            default:
                STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Values have unsupported type " << sourceType << ".");
        }
    }

    template<typename ValueType, typename SourceTypeEnum>
        requires std::is_enum_v<typename SourceTypeEnum::E>
    static std::vector<ValueType> createDecodedVector(storm::umb::GenericVector const& input, SourceTypeEnum sourceType, CSR const& csr = {}) {
        return applyDecodedVector<ValueType>(
            [](auto&& decodedInput) {
                std::vector<ValueType> v;
                v.reserve(std::ranges::size(decodedInput));
                for (auto&& value : decodedInput) {
                    v.push_back(value);
                }
                return v;  //{std::ranges::begin(decodedInput), std::ranges::end(decodedInput)};
            },
            input, sourceType, csr);
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, uint64_t>
    static auto uint64ToRationalRangeView(InputRange&& input) {
        STORM_LOG_ASSERT(std::ranges::size(input) % 2 == 0, "Input size is not even: " << std::ranges::size(input));
        return std::ranges::iota_view(0ull, std::ranges::size(input) / 2) | std::ranges::views::transform([&input](auto i) -> storm::RationalNumber {
                   // numerator is signed
                   return storm::utility::convertNumber<storm::RationalNumber>(static_cast<int64_t>(input[2 * i])) /
                          storm::utility::convertNumber<storm::RationalNumber>(input[2 * i + 1]);
               });
    }

    template<bool Signed, std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, uint64_t>
    static storm::NumberTraits<storm::RationalNumber>::IntegerType decodeArbitraryPrecisionInteger(InputRange&& input) {
        using IntegerType = typename storm::NumberTraits<storm::RationalNumber>::IntegerType;
        auto const twoTo64 = storm::utility::pow<IntegerType>(2, 64);

        STORM_LOG_ASSERT(std::ranges::size(input) > 0, "Input range must not be empty.");
        // We assume a little endian representation, so we reverse the input to start with the most significant bits.
        auto reverse_input = std::ranges::reverse_view(std::forward<InputRange>(input));

        // Helper function to decode from an unsigned range (with the most significant bits first).
        auto decodeFromUnsignedRange = [](auto&& input) -> IntegerType {
            auto const twoTo64 = storm::utility::pow<IntegerType>(2, 64);
            auto inputIt = std::ranges::begin(input);
            auto const inputEnd = std::ranges::end(input);
            auto result = storm::utility::convertNumber<IntegerType>(*inputIt);
            for (++inputIt; inputIt != inputEnd; ++inputIt) {
                result *= twoTo64;
                result += storm::utility::convertNumber<IntegerType>(*inputIt);
            }
            return result;
        };

        if constexpr (Signed) {
            // Find out if the number is negative, in which case we would compute the two's complement.
            uint64_t constexpr mostSignificantBitMask = 1ull << 63;
            if (*std::ranges::begin(reverse_input) & mostSignificantBitMask) {
                // Two's complement (e.g. 1111...1101 is -3)
                return -decodeFromUnsignedRange(reverse_input | std::ranges::views::transform([](uint64_t value) { return ~value; })) - 1;
            }
        }
        // Reaching this point means that the number is positive (i.e. signed and unsigned representations coincide).
        return decodeFromUnsignedRange(reverse_input);
    }

    template<std::ranges::input_range InputRange, std::ranges::input_range CsrRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, uint64_t> and std::same_as<std::ranges::range_value_t<CsrRange>, uint64_t>
    static auto uint64ToRationalRangeView(InputRange&& input, CsrRange&& csr) {
        STORM_LOG_ASSERT(std::ranges::size(csr) > 0, "CSR must not be empty.");
        auto const numEntries = std::ranges::size(csr) - 1;
        STORM_LOG_ASSERT(std::ranges::size(csr) % 2 == 0, "Input size is not even: " << std::ranges::size(input));
        return std::ranges::iota_view(0ull, numEntries) | std::ranges::views::transform([&input, csr](auto i) -> storm::RationalNumber {
                   auto const left = 2 * csr[i];
                   auto const right = 2 * csr[i + 1];
                   auto mid = left + (right - left) / 2;

                   std::ranges::subrange numeratorRange{std::ranges::begin(input) + left, std::ranges::begin(input) + mid};
                   storm::RationalNumber numerator = decodeArbitraryPrecisionInteger<true>(numeratorRange);  // signed

                   std::ranges::subrange denominatorRange{std::ranges::begin(input) + mid, std::ranges::begin(input) + right};
                   storm::RationalNumber denominator = decodeArbitraryPrecisionInteger<false>(denominatorRange);  // unsigned

                   return numerator / denominator;
               });
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, storm::RationalNumber>
    static bool rationalVectorRequiresCsr(InputRange&& input) {
        return std::any_of(std::ranges::begin(input), std::ranges::end(input), [](storm::RationalNumber const& r) {
            return storm::utility::numerator(r) < storm::utility::convertNumber<storm::RationalNumber>(std::numeric_limits<int64_t>::min()) ||
                   storm::utility::numerator(r) > storm::utility::convertNumber<storm::RationalNumber>(std::numeric_limits<int64_t>::max()) ||
                   storm::utility::denominator(r) > storm::utility::convertNumber<storm::RationalNumber>(std::numeric_limits<uint64_t>::max());
        });
    }
    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, storm::RationalNumber>
    static auto rationalToUint64ViewNoCsr(InputRange&& input) {
        static_assert(storm::RationalNumberDenominatorAlwaysPositive);
        return std::ranges::iota_view(0ull, std::ranges::size(input) * 2) | std::views::transform([&input](auto i) -> uint64_t {
                   storm::RationalNumber const& r = input[i / 2];
                   if (i % 2 == 0) {
                       // numerator
                       return static_cast<uint64_t>(storm::utility::convertNumber<int64_t, storm::RationalNumber>(storm::utility::numerator(r)));
                   } else {
                       // denominator
                       // We may assume that the denominator is always positive as this is a requirement for both GMP and CLN
                       STORM_LOG_ASSERT(storm::utility::denominator(r) > 0, "Denominator is not positive: " << storm::utility::denominator(r));
                       return storm::utility::convertNumber<uint64_t, storm::RationalNumber>(storm::utility::denominator(r));
                   }
               });
    }

    template<bool Signed>
    static uint64_t appendEncodedInteger(std::vector<uint64_t>& result, typename storm::NumberTraits<storm::RationalNumber>::IntegerType const& value) {
        if constexpr (Signed) {
            if (value < 0) {
                // Two's complement representation (e.g. -3 is 1111...1101)
                // We encode the corresponding positive value and then invert the bits.
                auto buckets = appendEncodedInteger<true>(result, -(value + 1));
                // Invert the bits
                for (auto& v : std::span<uint64_t>(result.end() - buckets, buckets)) {
                    v = ~v;
                }
                return buckets;
            } else {
                // We encode the non-negative value as if it were unsigned.
                auto buckets = appendEncodedInteger<false>(result, value);
                // Special case: the most significant bit must not be set. Otherwise, it would indicate a negative number in two's complement.
                uint64_t constexpr mostSignificantBitMask = 1ull << 63;
                if (result.back() & mostSignificantBitMask) {
                    result.push_back(0);
                    ++buckets;
                }
                return buckets;
            }
        } else {
            using IntegerType = typename storm::NumberTraits<storm::RationalNumber>::IntegerType;
            auto const twoTo64 = storm::utility::pow<IntegerType>(2, 64);
            // We assume a little endian representation, so we start with the least significant bits.

            STORM_LOG_ASSERT(value >= 0, "Value must be non-negative for unsigned encoding.");
            auto divisionResult = storm::utility::divide<IntegerType>(value, twoTo64);
            result.push_back(storm::utility::convertNumber<uint64_t, storm::RationalNumber>(divisionResult.second));
            uint64_t buckets = 1;
            while (divisionResult.first != 0) {
                divisionResult = storm::utility::divide<IntegerType>(divisionResult.first, twoTo64);
                result.push_back(storm::utility::convertNumber<uint64_t, storm::RationalNumber>(divisionResult.second));
                ++buckets;
            }
            return buckets;
        }
    }

    static uint64_t appendEncodedRational(std::vector<uint64_t>& result, storm::RationalNumber const& value) {
        // We may assume that the denominator is always positive as this is a requirement for both GMP and CLN
        static_assert(storm::RationalNumberDenominatorAlwaysPositive);
        auto const numeratorBuckets = appendEncodedInteger<true>(result, storm::utility::numerator(value));       // signed
        auto const denominatorBuckets = appendEncodedInteger<false>(result, storm::utility::denominator(value));  // unsigned

        if (numeratorBuckets < denominatorBuckets) {
            // add numerator buckets, requiring to move the denominator buckets to the end
            auto const oldStartOfDenominator = result.size() - denominatorBuckets;
            uint64_t const bucketsToAdd = denominatorBuckets - numeratorBuckets;
            auto const newStartOfDenominator = oldStartOfDenominator + bucketsToAdd;

            // move the denominator buckets to the new position
            result.resize(result.size() + bucketsToAdd);
            uint64_t i = result.size() - 1;
            for (; i >= newStartOfDenominator; --i) {
                result[i] = result[i - bucketsToAdd];
            }
            // fill the added numerator buckets with zeros
            for (; i >= oldStartOfDenominator; --i) {
                result[i] = 0;
            }
            return denominatorBuckets * 2;
        } else if (numeratorBuckets > denominatorBuckets) {
            // add denominator buckets to the end
            result.resize(result.size() + numeratorBuckets - denominatorBuckets, 0);
            return numeratorBuckets * 2;
        } else {
            // numeratorBuckets == denominatorBuckets
            return numeratorBuckets * 2;
        }
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, storm::RationalNumber>
    static std::pair<std::vector<uint64_t>, std::vector<uint64_t>> createUint64AndCsrFromRationalRange(InputRange&& input) {
        std::vector<uint64_t> values, csr;
        csr.reserve(std::ranges::size(input) + 1);
        csr.push_back(0);
        values.reserve(std::ranges::size(input) * 4);  // we guess that on average 256 bits are required per rational number...
        for (auto const& r : input) {
            appendEncodedRational(values, r);
            csr.push_back(values.size() / 2);  // we store the number of 128-bit buckets used for the numerator and denominator
        }
        values.shrink_to_fit();
        csr.shrink_to_fit();
        return {values, csr};
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, double>
    static auto doubleToIntervalRangeView(InputRange&& input) {
        STORM_LOG_ASSERT(std::ranges::size(input) % 2 == 0, "Input size is not even: " << std::ranges::size(input));
        return std::ranges::iota_view(0ull, std::ranges::size(input) / 2) |
               std::views::transform([&input](auto i) -> storm::Interval { return storm::Interval{input[2 * i], input[2 * i + 1]}; });
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, storm::Interval>
    static auto intervalToDoubleRangeView(InputRange&& input) {
        return std::ranges::iota_view(0ull, std::ranges::size(input) * 2) | std::views::transform([&input](auto i) {
                   if (i % 2 == 0) {
                       return storm::utility::convertNumber<double>(input[i / 2].lower());
                   } else {
                       return storm::utility::convertNumber<double>(input[i / 2].upper());
                   }
               });
    }

   private:
    template<typename TargetType, std::ranges::range Range>
    static auto conversionView(Range&& input) {
        using SourceType = std::ranges::range_value_t<Range>;
        if constexpr (std::same_as<TargetType, SourceType>) {
            return input;
        } else {
            return input | std::ranges::views::transform(
                               [](SourceType const& value) -> TargetType { return storm::utility::convertNumber<TargetType, SourceType>(value); });
        }
    }
};

}  // namespace storm::umb
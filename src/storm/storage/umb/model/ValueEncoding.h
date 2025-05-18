#pragma once
#include <ranges>
#include <span>
#include <vector>

#include "storm/exceptions/NotSupportedException.h"
#include "storm/storage/umb/model/GenericVector.h"
#include "storm/storage/umb/model/VectorType.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

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
    template<typename ValueType, StorageType Storage, typename SourceTypeEnum>
        requires std::is_enum_v<typename SourceTypeEnum::E>
    static auto applyDecodedVector(auto&& func, storm::umb::GenericVector<Storage> const& input, SourceTypeEnum sourceType,
                                   std::optional<VectorType<uint64_t, Storage>> const& csr = {}) {
        using E = decltype(sourceType)::E;
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

    template<typename ValueType, StorageType Storage, typename SourceTypeEnum>
        requires std::is_enum_v<typename SourceTypeEnum::E>
    static std::vector<ValueType> createDecodedVector(storm::umb::GenericVector<Storage> const& input, SourceTypeEnum sourceType,
                                                      std::optional<VectorType<uint64_t, Storage>> const& csr = {}) {
        return applyDecodedVector<ValueType, Storage>(
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

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, uint64_t>
    static auto uint64ToRationalRangeView(InputRange&& input, InputRange&& csr) {
        STORM_LOG_ASSERT(std::ranges::size(input) % 2 == 0, "Input size is not even: " << std::ranges::size(input));
        return std::ranges::iota_view(0ull, std::ranges::size(input) / 2) | std::ranges::views::transform([&input](auto i) -> storm::RationalNumber {
                   // numerator is signed
                   return storm::utility::convertNumber<storm::RationalNumber>(static_cast<int64_t>(input[2 * i])) /
                          storm::utility::convertNumber<storm::RationalNumber>(input[2 * i + 1]);
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
    static auto rationalToUint64View(InputRange&& input) {
        return std::ranges::iota_view(0ull, std::ranges::size(input) * 2) | std::views::transform([&input](auto i) -> uint64_t {
                   storm::RationalNumber const& r = input[i / 2];
                   if (i % 2 == 0) {
                       // numerator
                       return static_cast<uint64_t>(storm::utility::convertNumber<int64_t, storm::RationalNumber>(storm::utility::numerator(r)));
                   } else {
                       // denominator (we expect that they are always positive)
                       STORM_LOG_ASSERT(storm::utility::denominator(r) > 0, "Denominator is not positive: " << storm::utility::denominator(r));
                       return storm::utility::convertNumber<uint64_t, storm::RationalNumber>(storm::utility::denominator(r));
                   }
               });
    }

    template<std::ranges::input_range InputRange>
        requires std::same_as<std::ranges::range_value_t<InputRange>, double>
    static auto doubleToIntervalRangeView(InputRange&& input) {
        STORM_LOG_ASSERT(std::ranges::size(input) % 2 == 0, "Input size is not even: " << std::ranges::size(input));
        return std::ranges::iota_view(0ull, std::ranges::size(input) / 2) | std::views::transform([&input](auto i) -> storm::Interval {
                   return storm::Interval{input[2 * i], input[2 * i + 1]};
               });
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
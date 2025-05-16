#pragma once

#include <ranges>
#include <variant>

#include "storm/adapters/RationalNumberForward.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/VectorType.h"
#include "storm/utility/constants.h"

namespace storm::umb {

template<StorageType Storage>
class GenericVector {
   public:
    template<typename T>
    using Vec = VectorType<T, Storage>;

    template<typename T>
    void set(Vec<T>&& v) {
        data = std::move(v);
    }

    template<typename T>
    void set(Vec<T> const& v) {
        data = v;
    }

    void unset() {
        data = std::monostate();
    }

    template<typename T>
    Vec<T>& get() {
        return std::get<Vec<T>>(data);
    }

    template<typename T>
    Vec<T> const& get() const {
        return std::get<Vec<T>>(data);
    }

    template<typename T>
    bool isType() const {
        return std::holds_alternative<Vec<T>>(data);
    }

    bool hasValue() const {
        return !std::holds_alternative<std::monostate>(data);
    }

    template<typename FromType, typename ToType>
    auto convertFromTo() const {
        if constexpr (std::is_same_v<FromType, ToType>) {
            return get<ToType>();
        } else {
            return get<FromType>() |
                   std::ranges::views::transform([](FromType const& value) -> ToType { return storm::utility::convertNumber<ToType, FromType>(value); });
        }
    }

    template<typename T>
    std::vector<T> asVector() const {
        if (isType<T>()) {
            return std::vector<T>(get<T>().begin(), get<T>().end());
        } else {
            if constexpr (!std::is_same_v<T, bool>) {
                //                if (isType<int32_t>()) {
                //                    auto convertedView = convertFromTo<int32_t, T>();
                //                    return std::vector<T>(convertedView.begin(), convertedView.end());
                if (isType<uint64_t>()) {
                    auto convertedView = convertFromTo<uint64_t, T>();
                    return std::vector<T>(convertedView.begin(), convertedView.end());

                } else if (isType<double>()) {
                    auto convertedView = convertFromTo<double, T>();
                    return std::vector<T>(convertedView.begin(), convertedView.end());
                } else if (isType<storm::RationalNumber>()) {
                    auto convertedView = convertFromTo<storm::RationalNumber, T>();
                    return std::vector<T>(convertedView.begin(), convertedView.end());
                }
            }
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected type.");
        }
    }

    template<typename ValueType>
    ValueType at(uint64_t index) const {
        STORM_LOG_ASSERT(isType<double>() || isType<storm::RationalNumber>(), "unexpected type");
        if (isType<double>()) {
            return storm::utility::convertNumber<ValueType>(get<double>()[index]);
        } else {
            return storm::utility::convertNumber<ValueType>(get<storm::RationalNumber>()[index]);
        }
    }

   private:
    std::variant<std::monostate, Vec<bool>, Vec<int32_t>, Vec<uint64_t>, Vec<double>, Vec<storm::RationalNumber>, Vec<storm::Interval>> data;
};
}  // namespace storm::umb
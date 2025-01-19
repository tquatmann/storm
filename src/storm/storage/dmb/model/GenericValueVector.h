#include <ranges>
#include <variant>

#include "storm/adapters/RationalNumberForward.h"
#include "storm/storage/dmb/model/StorageType.h"
#include "storm/storage/dmb/model/VectorType.h"
#include "storm/utility/constants.h"

namespace storm::dmb {

template<StorageType Storage>
class GenericValueVectorType {
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

    template<typename FromType, typename ToType>
    auto convertFromTo() const {
        if constexpr (std::is_same_v<FromType, ToType>) {
            return get<ToType>();
        } else {
            return get<FromType>() | std::ranges::transform_view([](FromType const& value) { return storm::utility::convertNumber<ToType>(value); });
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
    std::variant<std::monostate, Vec<double>, Vec<storm::RationalNumber>> data;
};
}  // namespace storm::dmb
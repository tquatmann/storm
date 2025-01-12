#pragma once

#include "storm/adapters/JsonForward.h"

#include <boost/pfr.hpp>
#include "storm/exceptions/WrongFormatException.h"
#include "storm/utility/macros.h"

namespace storm {
/// Helper struct to enable json serialization
struct JsonSerialization {};
}  // namespace storm

NLOHMANN_JSON_NAMESPACE_BEGIN

template<typename T>
concept StormEnableSerializationConcept = std::is_same_v<typename T::JsonSerialization, ::storm::JsonSerialization>;

template<typename T>
concept StormIsOptional = std::same_as<T, std::optional<typename T::value_type>>;

template<StormEnableSerializationConcept T>
struct adl_serializer<T> {
    static_assert(T::JsonKeys.size() == boost::pfr::tuple_size_v<T>, "Number of JsonKeys does not match number of fields in struct.");

    template<typename JsonType>
    static void to_json(JsonType& json, T const& val) {
        to_json_i<0>(json, val);
        if (json.is_null()) {
            json = typename JsonType::object_t();
        }
    }

    template<typename JsonType>
    static void from_json(JsonType const& json, T& val) {
        STORM_LOG_THROW(json.is_object(), ::storm::exceptions::WrongFormatException, "Expected an object, got something else.");
        from_json_i<0>(json, val);
        // TODO: check if additional fields are present
    }

   private:
    static auto get_key(std::size_t i) {
        return std::data(T::JsonKeys)[i];
    }

    template<std::size_t I, typename JsonType>
    static void to_json_i(JsonType& json, T const& val) {
        if constexpr (I < boost::pfr::tuple_size_v<T>) {
            if constexpr (StormIsOptional<boost::pfr::tuple_element_t<I, T>>) {
                if (boost::pfr::get<I>(val).has_value()) {
                    json[get_key(I)] = boost::pfr::get<I>(val).value();
                }
            } else {
                json[get_key(I)] = boost::pfr::get<I>(val);
            }
            to_json_i<I + 1>(json, val);
        }
    }

    template<std::size_t I, typename JsonType>
    static void from_json_i(JsonType const& json, T& val) {
        if constexpr (I < boost::pfr::tuple_size_v<T>) {
            if constexpr (StormIsOptional<boost::pfr::tuple_element_t<I, T>>) {
                if (json.contains(get_key(I))) {
                    boost::pfr::get<I>(val) = json.at(get_key(I)).template get<typename boost::pfr::tuple_element_t<I, T>::value_type>();
                } else {
                    boost::pfr::get<I>(val) = std::nullopt;
                }
            } else {
                json.at(get_key(I)).get_to(boost::pfr::get<I>(val));
            }
            from_json_i<I + 1>(json, val);
        }
    }
};
NLOHMANN_JSON_NAMESPACE_END
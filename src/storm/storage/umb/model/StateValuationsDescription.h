#pragma once

#include <cstdint>
#include <optional>
#include <variant>

#include "storm/adapters/JsonAdapter.h"
#include "storm/adapters/JsonSerializationAdapter.h"

namespace storm::umb {

struct StateValuationsDescription {
    uint64_t alignment;
    struct Padding {
        uint64_t padding{0};
        static auto constexpr JsonKeys = {"padding"};
        using JsonSerialization = storm::JsonSerialization;
    };
    struct Variable {
        std::string name;
        std::optional<uint64_t> size;
        enum class Type { Bool, Uint, Int, Double, Rational, String };
        struct TypeDeclaration {
            using Values = Type;
            auto static constexpr Keys = {"bool", "uint", "int", "double", "rational", "string"};
        };
        storm::SerializedEnum<TypeDeclaration> type;
        std::optional<int64_t> lower, upper, offset;
        static auto constexpr JsonKeys = {"name", "size", "type", "lower", "upper", "offset"};
        using JsonSerialization = storm::JsonSerialization;
    };
    std::vector<std::variant<Padding, Variable>> variables;
    static auto constexpr JsonKeys = {"alignment", "variables"};
    using JsonSerialization = storm::JsonSerialization;
};
}  // namespace storm::umb

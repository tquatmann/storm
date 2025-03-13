#pragma once

#include <cstdint>
#include <ctime>
#include <optional>

#include "storm/adapters/JsonAdapter.h"
#include "storm/adapters/JsonSerializationAdapter.h"

namespace storm::umb {

struct ModelIndex {
    uint64_t formatVersion{0}, formatRevision{0};

    struct ModelData {
        std::optional<std::string> name, version;
        std::optional<std::vector<std::string>> authors;
        std::optional<std::string> description, comment, doi, url;
        static auto constexpr JsonKeys = {"name", "version", "authors", "description", "comment", "doi", "url"};
        using JsonSerialization = storm::JsonSerialization;
    };
    std::optional<ModelData> modelData;

    struct FileData {
        std::optional<std::string> tool, version;
        std::optional<uint64_t> date;
        std::optional<std::string> parameters;
        auto static constexpr JsonKeys = {"tool", "tool-version", "creation-date", "parameters"};
        using JsonSerialization = storm::JsonSerialization;

        std::string dateAsString() const {
            if (date) {
                std::time_t t = static_cast<std::time_t>(date.value());
                return std::ctime(&t);
            } else {
                return "unknown";
            }
        }
        void setDateToNow() {
            date = static_cast<uint64_t>(std::time(nullptr));
        }
    };
    std::optional<FileData> fileData;

    struct TransitionSystem {
        enum class Time { Discrete, Stochastic, UrgentStochastic };
        struct TimeDeclaration {
            using Values = Time;
            auto static constexpr Keys = {"discrete", "stochastic", "urgent-stochastic"};
        };
        storm::SerializedEnum<TimeDeclaration> time;

        uint64_t players{0};

        enum class BranchValues { None, Number, Interval };
        struct BranchValuesDeclaration {
            using Values = BranchValues;
            auto static constexpr Keys = {"none", "number", "interval"};
        };
        storm::SerializedEnum<BranchValuesDeclaration> branchValues;

        uint64_t numStates{0}, numInitialStates{0}, numChoices{0}, numActions{1} /* TODO */, numBranches{0};

        enum class BranchValueType { Double, Rational };
        struct BranchValueTypeDeclaration {
            using Values = BranchValueType;
            auto static constexpr Keys = {"double", "rational"};
        };
        std::optional<storm::SerializedEnum<BranchValueTypeDeclaration>> branchValueType;

        auto static constexpr JsonKeys = {"time",     "#players", "branch-values", "#states",          "#initial-states",
                                          "#choices", "#actions", "#branches",     "branch-value-type"};

        using JsonSerialization = storm::JsonSerialization;

    } transitionSystem;

    struct Annotation {
        std::optional<std::string> name;

        enum class AppliesTo { States, Choices, Branches };
        struct AppliesToDeclaration {
            using Values = AppliesTo;
            auto static constexpr Keys = {"states", "choices", "branches"};
        };
        storm::SerializedEnum<AppliesToDeclaration> appliesTo;

        enum class Type { Bool, Int32, Double, Rational, String };
        struct TypeDeclaration {
            using Values = Type;
            auto static constexpr Keys = {"bool", "int-32", "double", "rational", "string"};
        };
        storm::SerializedEnum<TypeDeclaration> type{Type::Bool};

        // TODO: #strings, lower, upper

        auto static constexpr JsonKeys = {"name", "applies-to", "type"};
        using JsonSerialization = storm::JsonSerialization;
    };
    std::map<std::string, Annotation> annotations;

    auto static constexpr JsonKeys = {"format-version", "format-revision", "model-data", "file-data", "transition-system", "annotations"};
    using JsonSerialization = storm::JsonSerialization;
};

}  // namespace storm::umb
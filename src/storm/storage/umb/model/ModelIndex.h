#pragma once

#include <cstdint>
#include <ctime>
#include <optional>

#include "storm/adapters/JsonAdapter.h"
#include "storm/adapters/JsonSerializationAdapter.h"

namespace storm::umb {

struct ModelIndex {
    uint64_t formatVersion{0}, formatRevision{0};  // TODO: set to 1

    struct ModelData {
        std::optional<std::string> name, version;
        std::optional<std::vector<std::string>> authors;
        std::optional<std::string> description, comment, doi, url;
        static auto constexpr JsonKeys = {"name", "version", "authors", "description", "comment", "doi", "url"};
        using JsonSerialization = storm::JsonSerialization;
    };
    std::optional<ModelData> modelData;

    struct FileData {
        std::optional<std::string> tool, toolVersion;
        std::optional<uint64_t> creationDate;
        std::optional<storm::json<storm::RationalNumber>> parameters;
        auto static constexpr JsonKeys = {"tool", "tool-version", "creation-date", "parameters"};
        using JsonSerialization = storm::JsonSerialization;

        std::string creationDateAsString() const {
            if (creationDate) {
                std::time_t t = static_cast<std::time_t>(creationDate.value());
                return std::ctime(&t);
            } else {
                return "unknown";
            }
        }
        void setCreationDateToNow() {
            if (std::time_t result; std::time(&result) == static_cast<std::time_t>(-1)) {
                STORM_LOG_WARN("No creation time is set for UMB index: unable to get the current time.");
                creationDate.reset();
            } else {
                creationDate.emplace(result);
            }
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

        auto static constexpr InvalidNumber = std::numeric_limits<uint64_t>::max();  // initial value for counts
        uint64_t numPlayers{InvalidNumber}, numStates{InvalidNumber}, numInitialStates{InvalidNumber}, numChoices{InvalidNumber}, numActions{InvalidNumber},
            numBranches{InvalidNumber};

        enum class BranchProbabilityType { None, Double, Rational, DoubleInterval, RationalInterval };
        struct BranchProbabilityTypeDeclaration {
            using Values = BranchProbabilityType;
            auto static constexpr Keys = {"none", "double", "rational", "double-interval", "rational-interval"};
        };
        storm::SerializedEnum<BranchProbabilityTypeDeclaration> branchProbabilityType;
        std::optional<storm::SerializedEnum<BranchProbabilityTypeDeclaration>> exitRateType;

        auto static constexpr JsonKeys = {
            "time", "#players", "#states", "#initial-states", "#choices", "#actions", "#branches", "branch-probability-type", "exit-rate-type"};
        using JsonSerialization = storm::JsonSerialization;

    } transitionSystem;

    struct Annotations {
        struct Annotation {
            static std::string getValidIdentifierFromAlias(std::string const& alias);

            std::optional<std::string> alias, description;
            std::optional<storm::RationalNumber> lower, upper;
            enum class AppliesTo { States, Choices, Branches };
            struct AppliesToDeclaration {
                using Values = AppliesTo;
                auto static constexpr Keys = {"states", "choices", "branches"};
            };
            std::optional<std::vector<storm::SerializedEnum<AppliesToDeclaration>>> appliesTo;
            enum class Type { Bool, Uint64, Uint64Interval, Int64, Int64Interval, Double, DoubleInterval, Rational, RationalInterval, String };
            struct TypeDeclaration {
                using Values = Type;
                auto static constexpr Keys = {
                    "bool", "uint64", "uint64-interval", "int64", "int64-interval", "double", "double-interval", "rational", "rational-interval", "string"};
            };
            std::optional<storm::SerializedEnum<TypeDeclaration>> type;
            auto static constexpr JsonKeys = {"alias", "description", "lower", "upper", "applies-to", "type"};
            using JsonSerialization = storm::JsonSerialization;
        };
        using AnnotationMap = std::map<std::string, Annotation>;
        std::optional<AnnotationMap> rewards, aps;

        std::optional<std::string> findRewardName(std::string const& id) const;
        std::optional<std::string> findAtomicPropositionName(std::string const& id) const;

        auto static constexpr JsonKeys = {"rewards", "aps"};
        using JsonSerialization = storm::JsonSerialization;
    } annotations;

    auto static constexpr JsonKeys = {"format-version", "format-revision", "model-data", "file-data", "transition-system", "annotations"};
    using JsonSerialization = storm::JsonSerialization;
};

}  // namespace storm::umb
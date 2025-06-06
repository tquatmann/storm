#pragma once

#include <cstdint>
#include <ctime>
#include <format>
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

        auto static constexpr JsonKeys = {"time", "#players", "#states", "#initial-states", "#choices", "#actions", "#branches", "branch-probability-type"};
        using JsonSerialization = storm::JsonSerialization;

    } transitionSystem;

    struct Annotations {
        enum class AppliesTo { States, Choices, Branches };
        struct AppliesToDeclaration {
            using Values = AppliesTo;
            auto static constexpr Keys = {"states", "choices", "branches"};
        };
        struct Reward {
            std::optional<std::string> alias, description;
            std::optional<storm::RationalNumber> lower, upper;
            enum class Type { Double, Rational, DoubleInterval, RationalInterval };
            struct TypeDeclaration {
                using Values = Type;
                auto static constexpr Keys = {"double", "rational", "double-interval", "rational-interval"};
            };
            storm::SerializedEnum<TypeDeclaration> type;
            std::optional<std::vector<storm::SerializedEnum<AppliesToDeclaration>>> appliesTo;  // todo check if this stays here
            auto static constexpr JsonKeys = {"alias", "description", "lower", "upper", "type", "applies-to"};
            using JsonSerialization = storm::JsonSerialization;
        };
        std::optional<std::map<std::string, Reward>> rewards;

        struct AtomicProposition {
            std::optional<std::string> alias, description;
            enum class Type { Bool };
            struct TypeDeclaration {
                using Values = Type;
                auto static constexpr Keys = {"bool"};
            };
            std::optional<storm::SerializedEnum<TypeDeclaration>> type;                         // todo: check if this stays here
            std::optional<std::vector<storm::SerializedEnum<AppliesToDeclaration>>> appliesTo;  // todo check if this stays here
            auto static constexpr JsonKeys = {"alias", "description", "type", "applies-to"};
            using JsonSerialization = storm::JsonSerialization;
        };
        std::optional<std::map<std::string, AtomicProposition>> atomicPropositions;

        static std::pair<std::string, std::optional<std::string>> getAllowedNameAndAlias(std::string const& inputName) {
            auto isAllowed = [](auto ch) { return (std::isalnum(ch) && !std::isupper(ch)) || ch == '_' || ch == '-'; };

            if (std::all_of(inputName.begin(), inputName.end(), isAllowed)) {
                return {inputName, std::nullopt};  // no alias needed
            } else {
                std::string newName;
                for (auto ch : inputName) {
                    if (isAllowed(ch)) {
                        newName += ch;
                    } else {
                        newName += "_0x" + std::format("{:x}", ch) + '_';
                    }
                }
                return {newName, inputName};  // alias is the original name
            }
        }

        template<typename MapType>
        static std::optional<std::string> findAnnotationName(std::optional<MapType> const& map, std::string const& id) {
            if (!map) {
                return {};
            }
            if (map->contains(id)) {
                return id;
            }
            for (auto const& [name, annotation] : map.value()) {
                if (annotation.alias == id) {
                    return name;
                }
            }
            return {};
        }

        std::optional<std::string> findRewardName(std::string const& id) const {
            return findAnnotationName(rewards, id);
        }

        std::optional<std::string> findAtomicPropositionName(std::string const& id) const {
            return findAnnotationName(atomicPropositions, id);
        }

        auto static constexpr JsonKeys = {"rewards", "aps"};
        using JsonSerialization = storm::JsonSerialization;
    } annotations;

    auto static constexpr JsonKeys = {"format-version", "format-revision", "model-data", "file-data", "transition-system", "annotations"};
    using JsonSerialization = storm::JsonSerialization;
};

}  // namespace storm::umb
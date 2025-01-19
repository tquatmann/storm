#pragma once

#include <cstdint>
#include <ctime>
#include <optional>

#include "storm/adapters/JsonAdapter.h"
#include "storm/adapters/JsonSerializationAdapter.h"

namespace storm::dmb {

struct ModelIndex {
    uint64_t formatVersion{0}, formatRevision{0};

    struct Metadata {
        std::optional<std::string> version;
        std::optional<std::vector<std::string>> authors;
        std::optional<std::string> description, comment, doi, url;
        static auto constexpr JsonKeys = {"version", "authors", "description", "comment", "doi", "url"};
        using JsonSerialization = storm::JsonSerialization;
    } metadata;

    struct Creation {
        std::string tool;
        std::optional<std::string> version;
        std::optional<uint64_t> date;
        std::optional<std::string> parameters;
        auto static constexpr JsonKeys = {"tool", "version", "date", "parameters"};
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

    } creation;

    enum class Time { Discrete, Stochastic, UrgentStochastic };
    NLOHMANN_JSON_SERIALIZE_ENUM(Time, {{Time::Discrete, "discrete"}, {Time::Stochastic, "stochastic"}, {Time::UrgentStochastic, "urgent-stochastic"}})
    Time time{Time::Discrete};
    uint64_t players{0};
    enum class BranchValues { None, Number, Interval };
    NLOHMANN_JSON_SERIALIZE_ENUM(BranchValues, {{BranchValues::None, "none"}, {BranchValues::Number, "number"}, {BranchValues::Interval, "interval"}})
    BranchValues branchValues{BranchValues::None};

    uint64_t numStates{0}, numInitialStates{0}, numChoices{0}, numBranches{0};

    enum class BranchValueType { Double, Rational };
    NLOHMANN_JSON_SERIALIZE_ENUM(BranchValueType, {{BranchValueType::Double, "double"}, {BranchValueType::Rational, "rational"}})
    BranchValueType branchValueType{BranchValueType::Double};

    auto static constexpr JsonKeys = {"format-version", "format-revision", "metadata",        "creation", "time",      "players",
                                      "branch-values",  "#states",         "#initial-states", "#choices", "#branches", "branch-value-type"};
    using JsonSerialization = storm::JsonSerialization;
};

}  // namespace storm::dmb
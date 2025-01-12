#pragma once

#include <cstdint>
#include <ctime>
#include <optional>

#include "storm/adapters/JsonSerializationAdapter.h"

namespace storm::dmb {

struct Index {
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

        std::string dateAsString() const {
            if (date) {
                std::time_t t = static_cast<std::time_t>(date.value());
                return std::ctime(&t);
            } else {
                return "unknown";
            }
        }

        auto static constexpr JsonKeys = {"tool", "version", "date", "parameters"};
        using JsonSerialization = storm::JsonSerialization;

    } creation;
    auto static constexpr JsonKeys = {"format-version", "format-revision", "metadata", "creation"};
    using JsonSerialization = storm::JsonSerialization;
};

}  // namespace storm::dmb
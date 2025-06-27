#include "storm/storage/umb/model/ModelIndex.h"

#include <algorithm>
#include <sstream>

namespace storm::umb {

std::string ModelIndex::Annotations::Annotation::getValidIdentifierFromAlias(std::string const& alias) {
    auto isAllowed = [](auto ch) { return (std::isalnum(ch) && !std::isupper(ch)) || ch == '_' || ch == '-'; };
    if (alias.empty()) {
        return "default";  // empty identifier is not allowed, so we return a default name
    }

    std::stringstream identifier;
    for (auto ch : alias) {
        if (isAllowed(ch)) {
            identifier << ch;
        } else if (ch == ' ') {  // replace spaces with underscores
            identifier << '_';
        } else {  // replace other characters with their hex representation
            identifier << "0x" << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(ch);
        }
    }
    return identifier.str();
}

template<typename MapType>
std::optional<std::string> findAnnotationName(std::optional<MapType> const& map, std::string const& id) {
    if (!map) {
        return {};
    }
    for (auto const& [name, annotation] : map.value()) {
        if (annotation.alias == id) {
            return name;
        }
    }
    if (map->contains(id)) {
        return id;
    }
    return {};
}

std::optional<std::string> ModelIndex::Annotations::findRewardName(std::string const& id) const {
    return findAnnotationName(rewards, id);
}

std::optional<std::string> ModelIndex::Annotations::findAtomicPropositionName(std::string const& id) const {
    return findAnnotationName(aps, id);
}
}  // namespace storm::umb
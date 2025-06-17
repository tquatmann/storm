#include "storm/storage/umb/model/ModelIndex.h"

#include <algorithm>
#include <sstream>

namespace storm::umb {

std::pair<std::string, std::optional<std::string>> ModelIndex::Annotations::getAllowedNameAndAlias(std::string const& inputName) {
    auto isAllowed = [](auto ch) { return (std::isalnum(ch) && !std::isupper(ch)) || ch == '_' || ch == '-'; };

    if (std::all_of(inputName.begin(), inputName.end(), isAllowed)) {
        return {inputName, inputName};  // no alias needed, but we still set it to the original name as a named annotation usually should have an alias
    } else {
        std::stringstream newName;
        for (auto ch : inputName) {
            if (isAllowed(ch)) {
                newName << ch;
            } else if (ch == ' ') {  // replace spaces with underscores
                newName << '_';
            } else {  // replace other characters with their hex representation
                newName << "0x" << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(ch);
            }
        }
        return {newName.str(), inputName};  // alias is the original name
    }
}
}  // namespace storm::umb
#include "AutomatonType.h"

#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

namespace storm::automata {
std::ostream& operator<<(std::ostream& stream, AutomatonType const& algorithm) {
    switch (algorithm) {
        case AutomatonType::DA:
            return stream << "da";
        case AutomatonType::LDBA:
            return stream << "ldba";
    }
    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unknown automaton type");
    return stream;
}

AutomatonType automatonTypeFromString(std::string const& type) {
    if (type == "da") {
        return AutomatonType::DA;
    } else if (type == "ldba") {
        return AutomatonType::LDBA;
    }
    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unknown automaton type: " << type);
}

}  // namespace storm::automata
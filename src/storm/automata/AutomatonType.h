#pragma once
#include <ostream>

namespace storm::automata {
enum class AutomatonType { DA, LDBA };

std::ostream& operator<<(std::ostream& stream, AutomatonType const& algorithm);
AutomatonType automatonTypeFromString(std::string const& algorithm);

}  // namespace storm::automata
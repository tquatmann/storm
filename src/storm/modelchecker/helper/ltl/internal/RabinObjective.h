#pragma once

#include <string>
#include <vector>

namespace storm::modelchecker::helper {

struct RabinCondition {
    std::vector<std::string> fin, inf;  // state labels that should be visited finitely / infinitely often (read as conjunction)
};
using RabinObjective = std::vector<RabinCondition>;  // disjunction of Rabin conditions
}  // namespace storm::modelchecker::helper
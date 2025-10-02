#pragma once

#include <string>
#include <vector>

namespace storm::modelchecker::helper {

struct RabinCondition {
    std::vector<std::string> inf;    // state labels that each should be visited infinitely often
    std::optional<std::string> fin;  // state label that should be visited only finitely often
};
using RabinObjective = std::vector<RabinCondition>;  // disjunction of Rabin conditions
}  // namespace storm::modelchecker::helper
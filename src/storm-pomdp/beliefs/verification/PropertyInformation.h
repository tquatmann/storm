#pragma once

#include <optional>
#include <set>
#include <string>

#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/solver/OptimizationDirection.h"

namespace storm::pomdp::beliefs {
struct PropertyInformation {
    enum class Kind { ReachabilityProbability, ExpectedTotalReachabilityReward };
    Kind kind;
    std::set<BeliefObservationType> targetObservations;
    std::optional<std::string> rewardModelName;
    storm::OptimizationDirection dir;
};
}  // namespace storm::pomdp::beliefs
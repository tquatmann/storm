#pragma once

#include <optional>

#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::bisimulation {

struct Options {
    enum class BisimulationType { Strong, Weak };

    BisimulationType bisimulationType = BisimulationType::Strong;

    std::optional<bool> preserveAllStateLabels;
    std::optional<bool> preserveAllRewards;
    bool preserveAllChoiceLabels = false;
    bool preserveChoiceOrigins = false;

    bool actionSensitive = false;

    std::optional<storm::RationalNumber> intervalApproximationEpsilon;
    std::optional<double> floatTolerance;
};

}  // namespace storm::bisimulation

#pragma once

#include <optional>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

struct Options {
    enum class BisimulationType { Strong, Weak };

    BisimulationType bisimulationType = BisimulationType::Strong;

    std::optional<bool> preserveAllStateLabels;
    std::optional<bool> preserveAllRewards;
    bool preserveAllChoiceLabels = false;
    bool preserveChoiceOrigins = false;

    bool actionSensitive = false;

    storm::RationalNumber tolerance = storm::utility::zero<RationalNumber>();
};

}  // namespace storm::bisimulation

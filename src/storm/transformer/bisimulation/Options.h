#pragma once

#include <optional>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

struct Options {
    enum class BisimulationType { Strong, Weak };

    // The model annotations that must be preserved.
    std::optional<bool> preserveAllStateLabels;  // if not specified, then all state labels are preserved iff no formula is given
    std::optional<bool> preserveAllRewards;      // if not specified, then all rewards are preserved iff no formula is given
    bool preserveChoiceLabels = false;
    bool preserveChoiceOrigins = false;

    // The kind of bisimulation that is applied.
    bool actionSensitive = false;  // if set, the i'th choice of state1 can only be matched with the i'th choice of state2
    BisimulationType bisimulationType = BisimulationType::Strong;
    storm::RationalNumber tolerance = storm::utility::zero<RationalNumber>();  // two values a and b are considered equal if |a-b| <= tolerance.

    // Algorithm Options
    bool createQuotientChoiceMapping = false;  // if set, a mapping from input choice index to quotient choice index is created and returned. This mapping can
                                               // be used, e.g. to map schedulers from the quotient model to the original model.
    bool preferSignatureRefinement =
        false;  // if set, signature-based refinement is used instead of splitter-based refinement (only applies to deterministic models)
};

}  // namespace storm::bisimulation

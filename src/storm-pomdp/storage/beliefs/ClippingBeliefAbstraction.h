#pragma once

#include <vector>

#include "storm-pomdp/storage/beliefs/BeliefBuilder.h"
#include "storm-pomdp/storage/beliefs/BeliefNumerics.h"
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"
#include "storm/solver/LpSolverForward.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
class ClippingBeliefAbstraction {
   public:
    using BeliefValueType = typename BeliefType::ValueType;

    struct BeliefClipping {
        bool isClippable;
        BeliefType targetBelief;
        BeliefValueType delta;
        BeliefFlatMap<BeliefValueType> deltaValues;
        bool onGrid = false;
    };

    ClippingBeliefAbstraction(std::vector<uint64_t>&& observationResolutions)
        : observationResolutions(std::forward<std::vector<uint64_t>>(observationResolutions)) {
        STORM_LOG_ASSERT(std::all_of(observationResolutions.begin(), observationResolutions.end(), [](auto o) { return o > 0; }),
                         "Expected that the resolutions are positive.");
    }

    template<typename AbstractCallback>
    void abstract(BeliefValueType&& probabilityFactor, BeliefType&& belief, AbstractCallback const& callback) {
        BeliefClipping clipping = clipBeliefToGrid(belief, observationResolutions[belief.observation()]);
        if (clipping.isClippable) {
            BeliefValueType a = (storm::utility::one<BeliefValueType>() - clipping.delta) * probabilityFactor;
            callback(std::move(a), std::move(clipping.targetBelief));
        } else {
            // Belief on Grid
            callback(std::move(probabilityFactor), std::move(belief));
        }
    }

    BeliefClipping clipBeliefToGrid(BeliefType const& belief, uint64_t resolution, const storm::storage::BitVector& isInfinite);

   private:
    std::vector<uint64_t> observationResolutions;
    std::shared_ptr<storm::solver::LpSolver<BeliefValueType>> lpSolver;
};

}  // namespace storm::pomdp::beliefs
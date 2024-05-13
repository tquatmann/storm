#pragma once

#include <vector>

#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/utility/BeliefNumerics.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/solver/LpSolverForward.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

/*!
 * @see 10.1007/978-3-030-99527-0_2
 */
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

    ClippingBeliefAbstraction(std::vector<uint64_t>&& observationResolutions);
    
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
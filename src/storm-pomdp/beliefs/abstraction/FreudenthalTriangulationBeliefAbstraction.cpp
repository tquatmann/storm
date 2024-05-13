#include "storm-pomdp/beliefs/abstraction/FreudenthalTriangulationBeliefAbstraction.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
FreudenthalTriangulationBeliefAbstraction<BeliefType>::FreudenthalTriangulationBeliefAbstraction(
    std::vector<BeliefValueType> const& observationTriangulationResolutions, FreudenthalTriangulationMode mode)
    : mode(mode) {
    observationResolutions.reserve(observationResolutions.size());
    for (auto const& res : observationTriangulationResolutions) {
        observationResolutions.emplace_back(storm::utility::ceil<BeliefValueType>(res));
        STORM_LOG_ASSERT(observationResolutions.back() > 0, "Expected that the resolution is a positive integer. Got " << res << " instead.");
    }
}

template class FreudenthalTriangulationBeliefAbstraction<Belief<double>>;
template class FreudenthalTriangulationBeliefAbstraction<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
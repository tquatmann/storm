#include "storm-pomdp/beliefs/abstraction/FreudenthalTriangulationBeliefAbstraction.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class FreudenthalTriangulationBeliefAbstraction<Belief<double>>;
template class FreudenthalTriangulationBeliefAbstraction<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
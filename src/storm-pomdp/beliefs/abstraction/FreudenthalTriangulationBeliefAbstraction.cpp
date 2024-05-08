#include "storm-pomdp/storage/beliefs/FreudenthalTriangulationBeliefAbstraction.h"
#include "storm-pomdp/storage/beliefs/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class FreudenthalTriangulationBeliefAbstraction<Belief<double>>;
template class FreudenthalTriangulationBeliefAbstraction<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
#include "storm-pomdp/storage/beliefs/BeliefBuilder.h"
#include "storm-pomdp/storage/beliefs/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class BeliefBuilder<Belief<double>>;
template class BeliefBuilder<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
#include "storm-pomdp/beliefs/storage/BeliefBuilder.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class BeliefBuilder<Belief<double>>;
template class BeliefBuilder<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
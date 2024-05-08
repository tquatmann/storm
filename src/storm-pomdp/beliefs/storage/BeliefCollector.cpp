#include "storm-pomdp/storage/beliefs/BeliefCollector.h"
#include "storm-pomdp/storage/beliefs/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class BeliefCollector<Belief<double>>;
template class BeliefCollector<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
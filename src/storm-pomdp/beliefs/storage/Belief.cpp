#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::pomdp::beliefs {

template class Belief<double>;
template class Belief<storm::RationalNumber>;

}  // namespace storm::pomdp::beliefs
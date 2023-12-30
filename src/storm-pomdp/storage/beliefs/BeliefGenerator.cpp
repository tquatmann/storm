#include "storm-pomdp/storage/beliefs/BeliefGenerator.h"
#include "storm-pomdp/storage/beliefs/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Pomdp.h"

namespace storm::pomdp::beliefs {

template class BeliefGenerator<storm::models::sparse::Pomdp<double>, Belief<double>>;
template class BeliefGenerator<storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<double>>;
template class BeliefGenerator<storm::models::sparse::Pomdp<double>, Belief<storm::RationalNumber>>;
template class BeliefGenerator<storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
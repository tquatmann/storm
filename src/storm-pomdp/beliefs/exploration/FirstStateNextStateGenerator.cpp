#include "storm-pomdp/beliefs/exploration/FirstStateNextStateGenerator.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Pomdp.h"

namespace storm::pomdp::beliefs {

template class FirstStateNextStateGenerator<storm::models::sparse::Pomdp<double>, Belief<double>>;
template class FirstStateNextStateGenerator<storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<double>>;
template class FirstStateNextStateGenerator<storm::models::sparse::Pomdp<double>, Belief<storm::RationalNumber>>;
template class FirstStateNextStateGenerator<storm::models::sparse::Pomdp<storm::RationalNumber>, Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
#pragma once

#include <functional>
#include <memory>

#include "storm-pomdp/beliefs/exploration/ExplorationInformation.h"
#include "storm-pomdp/beliefs/verification/PropertyInformation.h"
#include "storm/logic/Formulas.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"

namespace storm::pomdp::beliefs {

std::shared_ptr<storm::logic::Formula const> createFormulaForBeliefMdp(PropertyInformation const& propertyInformation);

// TODO: overloads for extra transition data (e.g. reward vectors)

template<typename BeliefMdpValueType, typename BeliefType>
std::shared_ptr<storm::models::sparse::Mdp<BeliefMdpValueType>> buildBeliefMdp(
    ExplorationInformation<BeliefMdpValueType, BeliefType> const& explorationInformation, PropertyInformation const& propertyInformation,
    std::function<BeliefMdpValueType(BeliefType const&)> computeCutOffValue);

}  // namespace storm::pomdp::beliefs
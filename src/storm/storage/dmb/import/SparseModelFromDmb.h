#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/dmb/model/DmbModelForward.h"

namespace storm::dmb {

storm::models::ModelType deriveModelType(storm::dmb::ModelIndex const& index);

template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> sparseModelFromDmb(storm::dmb::DmbModelBase const& dmbModel);

}  // namespace storm::dmb
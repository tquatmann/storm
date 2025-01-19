#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/storage/dmb/model/DmbModelForward.h"

namespace storm::dmb {
template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
std::unique_ptr<storm::dmb::DmbModelBase> sparseModelToDmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model);

}  // namespace storm::dmb
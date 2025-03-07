#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {
template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model);

}  // namespace storm::umb
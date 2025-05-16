#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {

storm::models::ModelType deriveModelType(storm::umb::ModelIndex const& index);

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModelFromUmb(storm::umb::UmbModelBase const& umbModel);

}  // namespace storm::umb
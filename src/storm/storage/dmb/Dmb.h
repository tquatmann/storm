#pragma once
#include <filesystem>
#include <memory>

#include "storm/storage/dmb/export/ExportOptions.h"
#include "storm/storage/dmb/import/ImportOptions.h"

#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"

namespace storm::dmb {
template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModelFromDmb(std::filesystem::path const& dmbLocation,
                                                                                            ImportOptions const& options = {});

template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
void exportModelToDmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model, std::filesystem::path const& targetLocation,
                      ExportOptions const& options = {});

}  // namespace storm::dmb

///*!
// * Load a model in dmb format from a file or directory and create the model.
// */
// template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
// std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModel(std::string const& file, DmbParserOptions const& options = {});
//}  // namespace storm::dmb
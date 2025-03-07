#pragma once
#include <filesystem>
#include <memory>

#include "storm/storage/umb/export/ExportOptions.h"
#include "storm/storage/umb/import/ImportOptions.h"

#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"

namespace storm::umb {
template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModelFromUmb(std::filesystem::path const& umbLocation,
                                                                                            ImportOptions const& options = {});

template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
void exportModelToUmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model, std::filesystem::path const& targetLocation,
                      ExportOptions const& options = {});

}  // namespace storm::umb

///*!
// * Load a model in umb format from a file or directory and create the model.
// */
// template<typename ValueType, typename RewardModelType = storm::models::sparse::StandardRewardModel<ValueType>>
// std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModel(std::string const& file, UmbParserOptions const& options = {});
//}  // namespace storm::umb
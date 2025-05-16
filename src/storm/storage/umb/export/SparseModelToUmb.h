#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"

#include "storm/storage/umb/export/ExportOptions.h"
#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {
template<typename ValueType>
storm::umb::UmbModelBase sparseModelToUmb(storm::models::sparse::Model<ValueType> const& model, ExportOptions const& options);

}  // namespace storm::umb
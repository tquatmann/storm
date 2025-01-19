#pragma once

#include <filesystem>
#include <memory>

#include "storm/storage/dmb/export/ExportOptions.h"
#include "storm/storage/dmb/model/DmbModelForward.h"

namespace storm::dmb {

void toDisk(storm::dmb::DmbModelBase const& dmbModel, std::filesystem::path const& dmbDir, ExportOptions const& options = {});

}  // namespace storm::dmb
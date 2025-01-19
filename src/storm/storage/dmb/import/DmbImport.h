#pragma once

#include <filesystem>
#include <memory>

#include "storm/storage/dmb/import/ImportOptions.h"
#include "storm/storage/dmb/model/DmbModelBase.h"

namespace storm::dmb {
class DmbModelBase;  // Forward Declaration

std::unique_ptr<storm::dmb::DmbModelBase> fromDisk(std::filesystem::path const& dmbDir, ImportOptions const& options = {});

}  // namespace storm::dmb
#pragma once

#include <filesystem>
#include <memory>

#include "storm/storage/umb/import/ImportOptions.h"
#include "storm/storage/umb/model/UmbModelBase.h"

namespace storm::umb {
class UmbModelBase;  // Forward Declaration

std::unique_ptr<storm::umb::UmbModelBase> fromDisk(std::filesystem::path const& umbDir, ImportOptions const& options = {});

}  // namespace storm::umb
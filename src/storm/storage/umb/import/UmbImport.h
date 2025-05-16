#pragma once

#include <filesystem>
#include <memory>

#include "storm/storage/umb/import/ImportOptions.h"
#include "storm/storage/umb/model/UmbModelBase.h"

namespace storm::umb {
class UmbModelBase;  // Forward Declaration

storm::umb::UmbModelBase importUmb(std::filesystem::path const& umbLocation, ImportOptions const& options = {});

}  // namespace storm::umb
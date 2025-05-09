#pragma once

#include <filesystem>
#include <memory>

#include "storm/storage/umb/export/ExportOptions.h"
#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {

void toDisk(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& options = {});
void toArchive(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& archivePath, ExportOptions const& options = {});

}  // namespace storm::umb
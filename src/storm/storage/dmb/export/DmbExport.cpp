#include "storm/storage/dmb/export/DmbExport.h"

#include "storm/storage/dmb/model/DmbModel.h"

#include "storm/exceptions/FileIoException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/io/BinaryFileWriter.h"
#include "storm/io/file.h"
#include "storm/utility/macros.h"

namespace storm::dmb {

namespace internal {

template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<typename VectorType::value_type>
void writeVectorToDisk(VectorType const& vector, std::filesystem::path const& dmbDir, std::filesystem::path const& fileName) {
    STORM_LOG_ASSERT(fileName.extension() == ".bin", "Undexpected file path '" << fileName << "'. File extension must be .bin");
    storm::io::BinaryFileWriter<typename VectorType::value_type, std::endian::little> writer(dmbDir / fileName);
    writer.write(vector);
}

template<typename VectorType>
void writeVectorToDisk(std::optional<VectorType> const& vector, std::filesystem::path const& dmbDir, std::filesystem::path const& fileName) {
    if (vector) {
        writeVectorToDisk(*vector, dmbDir, fileName);
    }
}

void writeIndexFile(storm::dmb::ModelIndex const& index, std::filesystem::path const& dmbDir) {
    std::ofstream stream;
    storm::json<storm::RationalNumber> indexJson(index);
    storm::io::openFile(dmbDir / "index.json", stream, false, true);
    stream << indexJson;
    storm::io::closeFile(stream);
}

template<StorageType Storage>
void writeStatesChoicesBranchesToDisk(DmbModel<Storage> const& dmbModel, std::filesystem::path const& dmbDir) {
    auto& states = dmbModel.states;
    writeVectorToDisk(states.stateToChoice, dmbDir, "state-to-choice.bin");
    writeVectorToDisk(states.stateToPlayer, dmbDir, "state-to-player.bin");
    // writeVectorToDisk(states.initialStates, dmbDir, "initial-states.bin");
    auto& choices = dmbModel.choices;
    writeVectorToDisk(choices.choiceToBranch, dmbDir, "choice-to-branch.bin");
    writeVectorToDisk(choices.choiceToAction, dmbDir, "choice-to-action.bin");
    // writeVectorToDisk(choices.actionStrings, dmbDir, "action-to-action-string.bin", "action-strings.bin");
    auto& branches = dmbModel.branches;
    writeVectorToDisk(branches.branchToTarget, dmbDir, "branch-to-target.bin");
    if (branches.branchToValue.template isType<double>()) {
        writeVectorToDisk(branches.branchToValue.template get<double>(), dmbDir, "branch-to-value.bin");
    } else if (branches.branchToValue.template isType<storm::RationalNumber>()) {
        // writeVectorToDisk(branches.branchToValue.template get<storm::RationalNumber>(), dmbDir, "branch-to-value.bin", "branch-rational.bin");
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected branch value type.");
    }
}

template<StorageType Storage>
void toDisk(storm::dmb::DmbModel<Storage> const& dmbModel, std::filesystem::path const& dmbDir, ExportOptions const& options) {
    std::filesystem::create_directories(dmbDir);
    writeIndexFile(dmbModel.index, dmbDir);
    writeStatesChoicesBranchesToDisk(dmbModel, dmbDir);
}

}  // namespace internal

void toDisk(storm::dmb::DmbModelBase const& dmbModel, std::filesystem::path const& dmbDir, ExportOptions const& options) {
    if (dmbModel.isStorageType(StorageType::Disk)) {
        internal::toDisk(dmbModel.template as<StorageType::Disk>(), dmbDir, options);
    } else if (dmbModel.isStorageType(StorageType::Memory)) {
        internal::toDisk(dmbModel.template as<StorageType::Memory>(), dmbDir, options);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

}  // namespace storm::dmb
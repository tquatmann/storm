#include "storm/storage/umb/export/UmbExport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/exceptions/FileIoException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/io/BinaryFileWriter.h"
#include "storm/io/file.h"
#include "storm/utility/macros.h"

namespace storm::umb {

namespace internal {

template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<typename VectorType::value_type>
void writeVectorToDisk(VectorType const& vector, std::filesystem::path const& umbDir, std::filesystem::path const& fileName) {
    STORM_LOG_ASSERT(fileName.extension() == ".bin", "Undexpected file path '" << fileName << "'. File extension must be .bin");
    storm::io::BinaryFileWriter<typename VectorType::value_type, std::endian::little> writer(umbDir / fileName);
    writer.write(vector);
}

void writeVectorToDisk(storm::storage::BitVector const& vector, std::filesystem::path const& umbDir, std::filesystem::path const& fileName) {
    STORM_LOG_ASSERT(fileName.extension() == ".bin", "Undexpected file path '" << fileName << "'. File extension must be .bin");
    storm::io::BinaryFileWriter<uint64_t, std::endian::little> writer(umbDir / fileName);
    for (uint64_t i = 0; i < vector.bucketCount(); ++i) {
        writer.write(vector.getBucket(i));
    }
}

template<StorageType Storage>
void writeVectorToDisk(UmbBitVector<Storage> const& vector, std::filesystem::path const& umbDir, std::filesystem::path const& fileName) {
    writeVectorToDisk(vector.getAsBitVectorAutoSize(), umbDir, fileName);
}

template<typename VectorType>
void writeVectorToDisk(std::optional<VectorType> const& vector, std::filesystem::path const& umbDir, std::filesystem::path const& fileName) {
    if (vector) {
        writeVectorToDisk(*vector, umbDir, fileName);
    }
}

void writeIndexFile(storm::umb::ModelIndex const& index, std::filesystem::path const& umbDir) {
    std::ofstream stream;
    storm::json<storm::RationalNumber> indexJson(index);
    storm::io::openFile(umbDir / "index.json", stream, false, true);
    stream << indexJson;
    storm::io::closeFile(stream);
}

template<StorageType Storage>
void writeStatesChoicesBranchesToDisk(UmbModel<Storage> const& umbModel, std::filesystem::path const& umbDir) {
    auto& states = umbModel.states;
    writeVectorToDisk(states.stateToChoice, umbDir, "state-to-choice.bin");
    writeVectorToDisk(states.stateToPlayer, umbDir, "state-to-player.bin");
    writeVectorToDisk(states.initialStates, umbDir, "initial-states.bin");
    auto& choices = umbModel.choices;
    writeVectorToDisk(choices.choiceToBranch, umbDir, "choice-to-branch.bin");
    writeVectorToDisk(choices.choiceToAction, umbDir, "choice-to-action.bin");
    // writeVectorToDisk(choices.actionStrings, umbDir, "action-to-action-string.bin", "action-strings.bin");
    auto& branches = umbModel.branches;
    writeVectorToDisk(branches.branchToTarget, umbDir, "branch-to-target.bin");
    if (branches.branchToValue.template isType<double>()) {
        writeVectorToDisk(branches.branchToValue.template get<double>(), umbDir, "branch-to-value.bin");
    } else if (branches.branchToValue.template isType<storm::RationalNumber>()) {
        // writeVectorToDisk(branches.branchToValue.template get<storm::RationalNumber>(), umbDir, "branch-to-value.bin", "branch-rational.bin");
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected branch value type.");
    }
}

template<StorageType Storage>
void toDisk(storm::umb::UmbModel<Storage> const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& options) {
    std::filesystem::create_directories(umbDir);
    writeIndexFile(umbModel.index, umbDir);
    writeStatesChoicesBranchesToDisk(umbModel, umbDir);
}

}  // namespace internal

void toDisk(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& options) {
    if (umbModel.isStorageType(StorageType::Disk)) {
        internal::toDisk(umbModel.template as<StorageType::Disk>(), umbDir, options);
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        internal::toDisk(umbModel.template as<StorageType::Memory>(), umbDir, options);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

}  // namespace storm::umb
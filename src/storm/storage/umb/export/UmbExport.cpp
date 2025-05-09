#include "storm/storage/umb/export/UmbExport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/io/ArchiveWriter.h"
#include "storm/io/BinaryFileWriter.h"
#include "storm/io/file.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::umb {

namespace detail {

template<typename T>
concept TargetType = std::same_as<std::remove_cvref_t<T>, std::filesystem::path> || std::same_as<std::remove_cvref_t<T>, storm::io::ArchiveWriter>;

void createDirectory(std::filesystem::path const& umbDir, std::filesystem::path const& subdirectory) {
    std::filesystem::create_directories(umbDir / subdirectory);
}

void createDirectory(storm::io::ArchiveWriter& archiveWriter, std::filesystem::path const& subdirectory) {
    archiveWriter.addDirectory(subdirectory.string());
}

/*!
 * Write a vector to disk.
 * The file path must have the extension .bin.
 */
template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<typename VectorType::value_type> || std::same_as<VectorType, storm::storage::BitVector> ||
             std::same_as<VectorType, storm::umb::UmbBitVector>
void writeVector(VectorType const& vector, std::filesystem::path const& umbDir, std::filesystem::path const& filepath) {
    STORM_LOG_ASSERT(filepath.extension() == ".bin", "Unexpected file path '" << filepath.filename() << "'. File extension must be .bin");
    if constexpr (std::is_same_v<VectorType, storm::umb::UmbBitVector>) {
        writeVector(vector.getAsBitVectorAutoSize(), umbDir, filepath);
    } else if constexpr (std::is_same_v<VectorType, storm::storage::BitVector>) {
        using BucketType = decltype(std::declval<storm::storage::BitVector&>().getBucket({}));
        storm::io::BinaryFileWriter<BucketType, std::endian::little> writer(umbDir / filepath);
        for (uint64_t i = 0; i < vector.bucketCount(); ++i) {
            writer.write(storm::utility::reverseBits(vector.getBucket(i)));
        }
    } else {
        storm::io::BinaryFileWriter<typename VectorType::value_type, std::endian::little> writer(umbDir / filepath);
        writer.write(vector);
    }
}

/*!
 * Write a vector to archive.
 * The file path must have the extension .bin.
 */
template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<typename VectorType::value_type> || std::same_as<VectorType, storm::storage::BitVector>
void writeVector(VectorType const& vector, storm::io::ArchiveWriter& archiveWriter, std::filesystem::path const& filepath) {
    STORM_LOG_ASSERT(filepath.extension() == ".bin", "Unexpected file path '" << filepath.filename() << "'. File extension must be .bin");
    archiveWriter.addBinaryFile(filepath.string(), vector);
}

template<StorageType Storage>
void writeVector(GenericVector<Storage> const& vector, TargetType auto& target, std::filesystem::path const& filepath) {
    if (vector.template isType<bool>()) {
        writeVector(vector.template get<bool>(), target, filepath);
    } else if (vector.template isType<int32_t>()) {
        writeVector(vector.template get<int32_t>(), target, filepath);
    } else if (vector.template isType<uint64_t>()) {
        writeVector(vector.template get<uint64_t>(), target, filepath);
    } else if (vector.template isType<double>()) {
        writeVector(vector.template get<double>(), target, filepath);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected type.");
    }
}

template<typename VectorType>
void writeVector(std::optional<VectorType> const& vector, TargetType auto& target, std::filesystem::path const& filepath) {
    if (vector) {
        writeVector(*vector, target, filepath);
    }
}

void writeIndexFile(storm::umb::ModelIndex const& index, std::filesystem::path const& umbDir) {
    std::ofstream stream;
    storm::json<storm::RationalNumber> indexJson(index);
    storm::io::openFile(umbDir / "index.json", stream, false, true);
    stream << storm::dumpJson(indexJson);
    storm::io::closeFile(stream);
}

void writeIndexFile(storm::umb::ModelIndex const& index, storm::io::ArchiveWriter& archiveWriter) {
    archiveWriter.addTextFile("index.json", storm::dumpJson(storm::json<storm::RationalNumber>(index)));
}

template<StorageType Storage>
void writeStatesChoicesBranches(UmbModel<Storage> const& umbModel, TargetType auto& target) {
    auto const& ts = umbModel.index.transitionSystem;
    auto const& states = umbModel.states;
    writeVector(states.stateToChoice, target, "state-to-choice.bin");
    writeVector(states.stateToPlayer, target, "state-to-player.bin");
    writeVector(states.initialStates, target, "initial-states.bin");
    auto const& choices = umbModel.choices;
    writeVector(choices.choiceToBranch, target, "choice-to-branch.bin");
    writeVector(choices.choiceToAction, target, "choice-to-action.bin");
    // writeVector(choices.actionStrings, target, "action-to-action-string.bin", "action-strings.bin");
    auto const& branches = umbModel.branches;
    writeVector(branches.branchToTarget, target, "branch-to-target.bin");

    using BranchValues = storm::umb::ModelIndex::TransitionSystem::BranchValues;
    using BranchValueType = storm::umb::ModelIndex::TransitionSystem::BranchValueType;
    if (ts.branchValues == BranchValues::Number) {
        STORM_LOG_THROW(ts.branchValueType.has_value(), storm::exceptions::WrongFormatException, "Branch values are numbers, but no type is specified.");
        if (auto const branchVT = *ts.branchValueType; branchVT == BranchValueType::Double) {
            STORM_LOG_ASSERT(branches.branchValues.template isType<double>(), "Unexpected branch value type.");
            writeVector(branches.branchValues, target, "branch-values.bin");
        } else {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value type " << branchVT << " unhandled.");
        }
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value kind " << ts.branchValues << " unhandled.");
    }
}

template<StorageType Storage>
void writeAnnotations(UmbModel<Storage> const& umbModel, TargetType auto& target) {
    createDirectory(target, "annotations");
    for (auto const& [annotationId, annotation] : umbModel.index.annotations) {
        std::filesystem::path const dirOfAnnotation = "annotations/" + annotationId;
        createDirectory(target, dirOfAnnotation);
        auto const& annotationFile = umbModel.annotations.at(annotationId);
        writeVector(annotationFile.values, target, dirOfAnnotation / "values.bin");
    }
}

template<StorageType Storage>
void exportUmb(storm::umb::UmbModel<Storage> const& umbModel, TargetType auto& target, ExportOptions const& options) {
    writeIndexFile(umbModel.index, target);
    writeStatesChoicesBranches(umbModel, target);
    writeAnnotations(umbModel, target);
}

}  // namespace detail

void toDisk(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& options) {
    std::filesystem::create_directories(umbDir);
    if (umbModel.isStorageType(StorageType::Disk)) {
        detail::exportUmb(umbModel.template as<StorageType::Disk>(), umbDir, options);
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        detail::exportUmb(umbModel.template as<StorageType::Memory>(), umbDir, options);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

void toArchive(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& archivePath, ExportOptions const& options) {
    storm::io::ArchiveWriter archiveWriter(archivePath);
    if (umbModel.isStorageType(StorageType::Disk)) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Exporting to archive from disk storage is not supported.");  // TODO
        //        detail::exportUmb(umbModel.template as<StorageType::Disk>(), archiveWriter, options);
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        detail::exportUmb(umbModel.template as<StorageType::Memory>(), archiveWriter, options);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

}  // namespace storm::umb
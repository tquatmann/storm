#include "storm/storage/umb/export/UmbExport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/exceptions/FileIoException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/io/BinaryFileWriter.h"
#include "storm/io/file.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::umb {

namespace internal {

template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<typename VectorType::value_type>
void writeVectorToDisk(VectorType const& vector, std::filesystem::path const& filepath) {
    STORM_LOG_ASSERT(filepath.extension() == ".bin", "Undexpected file path '" << filepath.filename() << "'. File extension must be .bin");
    storm::io::BinaryFileWriter<typename VectorType::value_type, std::endian::little> writer(filepath);
    writer.write(vector);
}

void writeVectorToDisk(storm::storage::BitVector const& vector, std::filesystem::path const& filepath) {
    STORM_LOG_ASSERT(filepath.extension() == ".bin", "Undexpected file path '" << filepath.filename() << "'. File extension must be .bin");
    storm::io::BinaryFileWriter<uint64_t, std::endian::little> writer(filepath);
    for (uint64_t i = 0; i < vector.bucketCount(); ++i) {
        writer.write(storm::utility::reverseBits(vector.getBucket(i)));
    }
}

template<StorageType Storage>
void writeVectorToDisk(UmbBitVector<Storage> const& vector, std::filesystem::path const& filepath) {
    writeVectorToDisk(vector.getAsBitVectorAutoSize(), filepath);
}

template<StorageType Storage>
void writeVectorToDisk(GenericVector<Storage> const& vector, std::filesystem::path const& filepath) {
    if (vector.template isType<bool>()) {
        writeVectorToDisk(vector.template get<bool>(), filepath);
    } else if (vector.template isType<int32_t>()) {
        writeVectorToDisk(vector.template get<int32_t>(), filepath);
    } else if (vector.template isType<uint64_t>()) {
        writeVectorToDisk(vector.template get<uint64_t>(), filepath);
    } else if (vector.template isType<double>()) {
        writeVectorToDisk(vector.template get<double>(), filepath);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected type.");
    }
}

template<typename VectorType>
void writeVectorToDisk(std::optional<VectorType> const& vector, std::filesystem::path const& filepath) {
    if (vector) {
        writeVectorToDisk(*vector, filepath);
    }
}

void writeIndexFile(storm::umb::ModelIndex const& index, std::filesystem::path const& umbDir) {
    std::ofstream stream;
    storm::json<storm::RationalNumber> indexJson(index);
    storm::io::openFile(umbDir / "index.json", stream, false, true);
    stream << storm::dumpJson(indexJson);
    storm::io::closeFile(stream);
}

template<StorageType Storage>
void writeStatesChoicesBranchesToDisk(UmbModel<Storage> const& umbModel, std::filesystem::path const& umbDir) {
    auto const& ts = umbModel.index.transitionSystem;
    auto const& states = umbModel.states;
    writeVectorToDisk(states.stateToChoice, umbDir / "state-to-choice.bin");
    writeVectorToDisk(states.stateToPlayer, umbDir / "state-to-player.bin");
    writeVectorToDisk(states.initialStates, umbDir / "initial-states.bin");
    auto const& choices = umbModel.choices;
    writeVectorToDisk(choices.choiceToBranch, umbDir / "choice-to-branch.bin");
    writeVectorToDisk(choices.choiceToAction, umbDir / "choice-to-action.bin");
    // writeVectorToDisk(choices.actionStrings, umbDir, "action-to-action-string.bin", "action-strings.bin");
    auto const& branches = umbModel.branches;
    writeVectorToDisk(branches.branchToTarget, umbDir / "branch-to-target.bin");

    using BranchValues = storm::umb::ModelIndex::TransitionSystem::BranchValues;
    using BranchValueType = storm::umb::ModelIndex::TransitionSystem::BranchValueType;
    if (ts.branchValues == BranchValues::Number) {
        STORM_LOG_THROW(ts.branchValueType.has_value(), storm::exceptions::WrongFormatException, "Branch values are numbers, but no type is specified.");
        if (auto const branchVT = *ts.branchValueType; branchVT == BranchValueType::Double) {
            STORM_LOG_ASSERT(branches.branchValues.template isType<double>(), "Unexpected branch value type.");
            writeVectorToDisk(branches.branchValues, umbDir / "branch-values.bin");
        } else {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value type " << branchVT << " unhandled.");
        }
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value kind " << ts.branchValues << " unhandled.");
    }
}

template<StorageType Storage>
void writeAnnotationsToDisk(UmbModel<Storage> const& umbModel, std::filesystem::path const& annotationsDir) {
    for (auto const& [annotationId, annotation] : umbModel.index.annotations) {
        auto const dirOfAnnotation = annotationsDir / annotationId;
        std::filesystem::create_directories(dirOfAnnotation);
        auto const& annotationFile = umbModel.annotations.at(annotationId);
        writeVectorToDisk(annotationFile.values, dirOfAnnotation / "values.bin");
    }
}

template<StorageType Storage>
void toDisk(storm::umb::UmbModel<Storage> const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& options) {
    std::filesystem::create_directories(umbDir);
    writeIndexFile(umbModel.index, umbDir);
    writeStatesChoicesBranchesToDisk(umbModel, umbDir);
    writeAnnotationsToDisk(umbModel, umbDir / "annotations");
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
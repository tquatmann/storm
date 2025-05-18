#include "storm/storage/umb/export/UmbExport.h"

#include <boost/pfr.hpp>

#include "storm/storage/umb/model/UmbModel.h"
#include "storm/storage/umb/model/ValueEncoding.h"

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
    requires storm::io::IsBinaryFileWritable<std::ranges::range_value_t<VectorType>> || std::same_as<VectorType, storm::storage::BitVector> ||
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
        storm::io::BinaryFileWriter<std::ranges::range_value_t<VectorType>, std::endian::little> writer(umbDir / filepath);
        writer.write(vector);
    }
}

/*!
 * Write a vector to archive.
 * The file path must have the extension .bin.
 */
template<typename VectorType>
    requires storm::io::IsBinaryFileWritable<std::ranges::range_value_t<VectorType>> || std::same_as<VectorType, storm::storage::BitVector>
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

void writeIndexFile(storm::umb::ModelIndex const& index, std::filesystem::path const& umbDir, std::filesystem::path const& filepath) {
    std::ofstream stream;
    storm::json<storm::RationalNumber> indexJson(index);
    storm::io::openFile(umbDir / filepath, stream, false, true);
    stream << storm::dumpJson(indexJson);
    storm::io::closeFile(stream);
}

void writeIndexFile(storm::umb::ModelIndex const& index, storm::io::ArchiveWriter& archiveWriter, std::filesystem::path const& filepath) {
    archiveWriter.addTextFile(filepath, storm::dumpJson(storm::json<storm::RationalNumber>(index)));
}

template<typename T>
concept HasFileNames = requires { T::FileNames.size(); };

template<typename T>
concept FileNameMap = std::same_as<std::remove_cvref_t<typename T::key_type>, std::string> && HasFileNames<typename T::mapped_type>;

template<typename T>
concept IsOptionalWithFileNames = std::same_as<std::remove_cvref_t<T>, std::optional<typename T::value_type>> && HasFileNames<typename T::value_type>;

template<typename UmbStructure>
void exportFiles(UmbStructure const& umbStructure, TargetType auto& target, std::filesystem::path const& context) {
    static_assert(UmbStructure::FileNames.size() == boost::pfr::tuple_size_v<UmbStructure>, "Number of file names does not match number of fields in struct.");
    boost::pfr::for_each_field(umbStructure, [&](auto const& field, std::size_t i) {
        // potentially create directory for sub-field
        std::filesystem::path fieldName = std::data(UmbStructure::FileNames)[i];
        // export the field, either with a recursive call or via writeVector
        using FieldType = std::remove_cvref_t<decltype(field)>;
        if constexpr (HasFileNames<FieldType>) {
            if (!fieldName.empty()) {
                createDirectory(target, context / fieldName);
            }
            exportFiles(field, target, context / fieldName);
        } else if constexpr (IsOptionalWithFileNames<FieldType>) {
            if (field) {
                if (!fieldName.empty()) {
                    createDirectory(target, context / fieldName);
                }
                exportFiles(*field, target, context / fieldName);
            }
        } else if constexpr (FileNameMap<FieldType>) {
            if (!field.empty() && !fieldName.empty()) {
                createDirectory(target, context / fieldName);
            }
            for (auto const& [key, value] : field) {
                createDirectory(target, context / fieldName / key);
                exportFiles(value, target, context / fieldName / key);
            }
        } else if constexpr (std::is_same_v<FieldType, storm::umb::ModelIndex>) {
            writeIndexFile(field, target, context / fieldName);
        } else if constexpr (std::is_same_v<FieldType, GenericVector<StorageType::Memory>> || std::is_same_v<FieldType, GenericVector<StorageType::Disk>>) {
            if (field.template isType<storm::RationalNumber>()) {
                if (ValueEncoding::rationalVectorRequiresCsr(field.template get<storm::RationalNumber>())) {
                    assert(false);  // todo
                } else {
                    writeVector(ValueEncoding::rationalToUint64View(field.template get<storm::RationalNumber>()), target, context / fieldName);
                }
            } else if (field.template isType<storm::Interval>()) {
                writeVector(ValueEncoding::intervalToDoubleRangeView(field.template get<storm::Interval>()), target, context / fieldName);
            } else {
                writeVector(field, target, context / fieldName);
            }
        } else {
            writeVector(field, target, context / fieldName);
        }
    });
}

}  // namespace detail

void toDisk(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& umbDir, ExportOptions const& /*options*/) {
    std::filesystem::create_directories(umbDir);
    if (umbModel.isStorageType(StorageType::Disk)) {
        detail::exportFiles(umbModel.template as<StorageType::Disk>(), umbDir, {});
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        detail::exportFiles(umbModel.template as<StorageType::Memory>(), umbDir, {});
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

void toArchive(storm::umb::UmbModelBase const& umbModel, std::filesystem::path const& archivePath, ExportOptions const& /*options*/) {
    storm::io::ArchiveWriter archiveWriter(archivePath);
    if (umbModel.isStorageType(StorageType::Disk)) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Exporting to archive from disk storage is not supported.");  // TODO
        //        detail::exportUmb(umbModel.template as<StorageType::Disk>(), archiveWriter);
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        detail::exportFiles(umbModel.template as<StorageType::Memory>(), archiveWriter, {});
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

}  // namespace storm::umb
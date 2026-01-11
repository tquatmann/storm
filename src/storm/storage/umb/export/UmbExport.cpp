#include "storm/storage/umb/export/UmbExport.h"

#include <boost/pfr.hpp>

#include "storm/storage/umb/model/UmbModel.h"
#include "storm/storage/umb/model/ValueEncoding.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/io/ArchiveWriter.h"
#include "storm/io/file.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::umb {

namespace detail {

template<typename T>

void createDirectory(std::filesystem::path const& umbDir, std::filesystem::path const& subdirectory) {
    std::filesystem::create_directories(umbDir / subdirectory);
}

void createDirectory(storm::io::ArchiveWriter& archiveWriter, std::filesystem::path const& subdirectory) {
    archiveWriter.addDirectory(subdirectory.string());
}

/*!
 * Write a vector to archive.
 * The file path must have the extension .bin.
 */
template<typename VectorType>
    requires(!std::is_same_v<std::remove_cvref_t<VectorType>, storm::umb::GenericVector>)
void writeVector(VectorType const& vector, storm::io::ArchiveWriter& archiveWriter, std::filesystem::path const& filepath) {
    STORM_LOG_ASSERT(filepath.extension() == ".bin", "Unexpected file path '" << filepath.filename() << "'. File extension must be .bin");
    archiveWriter.addBinaryFile(filepath.string(), vector);
}

void writeVector(GenericVector const& vector, storm::io::ArchiveWriter& target, std::filesystem::path const& filepath) {
    if (!vector.hasValue()) {
        return;
    }
    if (vector.template isType<bool>()) {
        writeVector(vector.template get<bool>(), target, filepath);
    } else if (vector.template isType<uint64_t>()) {
        writeVector(vector.template get<uint64_t>(), target, filepath);
    } else if (vector.template isType<int64_t>()) {
        writeVector(vector.template get<int64_t>(), target, filepath);
    } else if (vector.template isType<double>()) {
        writeVector(vector.template get<double>(), target, filepath);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected type.");
    }
}

template<typename VectorType>
void writeVector(std::optional<VectorType> const& vector, storm::io::ArchiveWriter& target, std::filesystem::path const& filepath) {
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
    requires HasFileNames<UmbStructure>
void applyToFieldWithName(UmbStructure const& umbStructure, auto&& fieldName, auto&& func) {
    auto constexpr i = std::find(UmbStructure::FileNames.begin(), UmbStructure::FileNames.end(), fieldName) - UmbStructure::FileNames.begin();
    static_assert(i < UmbStructure::FileNames.size(), "Field name not found in file names.");
    func(boost::pfr::get<i>(umbStructure));
}

template<typename UmbStructure>
    requires HasFileNames<UmbStructure>
void exportFiles(UmbStructure const& umbStructure, storm::io::ArchiveWriter& target, std::filesystem::path const& context) {
    static_assert(UmbStructure::FileNames.size() == boost::pfr::tuple_size_v<UmbStructure>, "Number of file names does not match number of fields in struct.");
    boost::pfr::for_each_field(umbStructure, [&](auto const& field, std::size_t fieldIndex) {
        // potentially create directory for sub-field
        std::filesystem::path fieldName = std::data(UmbStructure::FileNames)[fieldIndex];
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
        } else if constexpr (std::is_same_v<FieldType, GenericVector>) {
            if (field.template isType<storm::RationalNumber>()) {
                if (ValueEncoding::rationalVectorRequiresCsr(field.template get<storm::RationalNumber>())) {
                    // TODO: This is not very memory efficient as we create the entire vector in memory first instead of writing it in chunks.
                    auto [values, csr] = ValueEncoding::createUint64AndCsrFromRationalRange(field.template get<storm::RationalNumber>());
                    writeVector(std::move(values), target, context / fieldName);
                    std::string csrFieldName;
                    if (fieldName == "branch-probabilities.bin") {
                        csrFieldName = "branch-to-probability.bin";
                    } else if (fieldName == "exit-rates.bin") {
                        csrFieldName = "state-to-exit-rate.bin";
                    } else if (fieldName == "values.bin") {
                        csrFieldName = "to-value.bin";
                    } else {
                        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException,
                                        "Unexpected field name '" << fieldName << "' has no associated CSR file name.");
                    }
                    writeVector(std::move(csr), target, context / csrFieldName);
                } else {
                    writeVector(ValueEncoding::rationalToUint64ViewNoCsr(field.template get<storm::RationalNumber>()), target, context / fieldName);
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

void toArchive(storm::umb::UmbModel const& umbModel, std::filesystem::path const& archivePath, ExportOptions const& /*options*/) {
    storm::io::ArchiveWriter archiveWriter(archivePath);
    detail::exportFiles(umbModel, archiveWriter, {});
}

}  // namespace storm::umb
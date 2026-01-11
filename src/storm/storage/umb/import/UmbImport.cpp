#include "storm/storage/umb/import/UmbImport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/adapters/JsonAdapter.h"
#include "storm/io/ArchiveReader.h"
#include "storm/io/file.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/FileIoException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/WrongFormatException.h"

namespace storm::umb {

namespace internal {
template<typename T>
concept HasFileNames = requires { T::FileNames.size(); };

template<typename T>
concept FileNameMap = std::same_as<std::remove_cvref_t<typename T::key_type>, std::string> && HasFileNames<typename T::mapped_type>;

template<typename T>
concept IsOptional = std::same_as<std::remove_cvref_t<T>, std::optional<typename T::value_type>>;

template<typename T>
concept IsOptionalWithFileNames = IsOptional<T> && HasFileNames<typename T::value_type>;

void parseIndexFromDisk(std::filesystem::path const& indexFilePath, storm::umb::ModelIndex& index) {
    storm::json<storm::RationalNumber> parsedStructure;
    std::ifstream file;
    storm::io::openFile(indexFilePath, file);
    file >> parsedStructure;
    storm::io::closeFile(file);
    parsedStructure.get_to(index);
}

void parseIndexFromString(std::string const& indexFileString, storm::umb::ModelIndex& index) {
    storm::json<storm::RationalNumber>::parse(indexFileString).get_to(index);
}

/*!
 * Prepares annotations so that all fields are available according to the index.
 */
void prepareAnnotations(storm::umb::UmbModel& umbModel) {
    if (umbModel.index.annotations.rewards) {
        for (auto const& [name, rew] : umbModel.index.annotations.rewards.value()) {
            umbModel.rewards[name];
        }
    }
    if (umbModel.index.annotations.aps) {
        for (auto const& [name, ap] : umbModel.index.annotations.aps.value()) {
            umbModel.aps[name];
        }
    }
}

std::filesystem::path getFilePath(storm::io::ArchiveReadEntry const& src) {
    return src.name();
}

template<typename VecT>
    requires std::same_as<VecT, std::vector<typename VecT::value_type>>
VecT importVector(storm::io::ArchiveReadEntry& src) {
    return src.template toVector<typename VecT::value_type>();
}

template<typename VecT>
    requires std::same_as<VecT, storm::storage::BitVector>
VecT importVector(storm::io::ArchiveReadEntry& src) {
    return src.template toVector<bool>();
}

template<typename VecT>
void importVector(storm::io::ArchiveReadEntry& src, std::optional<VecT>& target) {
    target = importVector<VecT>(src);
}

template<typename VecT>
    requires(!IsOptional<VecT>)
void importVector(storm::io::ArchiveReadEntry& src, VecT& target) {
    target = importVector<VecT>(src);
}

template<typename ValueType>
void importGenericVector(storm::io::ArchiveReadEntry& src, GenericVector& target) {
    target.template set<ValueType>(importVector<typename GenericVector::template Vec<ValueType>>(src));
}

template<typename TypeDecl>
    requires std::is_enum_v<typename TypeDecl::E>
void importGenericVector(storm::io::ArchiveReadEntry& src, TypeDecl const& type, GenericVector& target) {
    using E = typename TypeDecl::E;
    // handle common types
    if (type == E::Double || type == E::DoubleInterval) {
        importGenericVector<double>(src, target);
    } else if (type == E::Rational || type == E::RationalInterval) {
        importGenericVector<uint64_t>(src, target);
    } else {
        if constexpr (std::is_same_v<E, storm::umb::ModelIndex::Annotations::Annotation::Type>) {
            // Handle annotation-specific types.
            if (type == storm::umb::ModelIndex::Annotations::Annotation::Type::Bool) {
                importGenericVector<bool>(src, target);
            } else if (type == storm::umb::ModelIndex::Annotations::Annotation::Type::Uint64) {
                importGenericVector<uint64_t>(src, target);
            } else if (type == storm::umb::ModelIndex::Annotations::Annotation::Type::Int64) {
                importGenericVector<int64_t>(src, target);
            } else {
                STORM_LOG_THROW(false, storm::exceptions::WrongFormatException,
                                "Annotation type " << type << " for vector located in '" << getFilePath(src) << "' is not handled");
            }
        } else {
            STORM_LOG_THROW(false, storm::exceptions::WrongFormatException,
                            "Type " << type << " for numeric vector located in '" << getFilePath(src) << "' is not handled");
        }
    }
}

void importGenericVector(storm::io::ArchiveReadEntry& src, storm::umb::ModelIndex const& index, GenericVector& target) {
    // Find type information in the index that matches the given src.
    auto srcPath = getFilePath(src);
    if (srcPath == "branch-probabilities.bin") {
        importGenericVector(src, index.transitionSystem.branchProbabilityType, target);
    } else if (srcPath == "exit-rates.bin") {
        auto const& rateType = index.transitionSystem.exitRateType;
        STORM_LOG_THROW(rateType.has_value(), storm::exceptions::WrongFormatException,
                        "Exit rate type is not set in the index, but exit rates are present in '" << srcPath << "'.");
        importGenericVector(src, rateType.value(), target);
    } else {
        // Reaching this point means that we must have annotation values.
        std::vector<std::string> p(srcPath.begin(), srcPath.end());
        // We expect a path of the form "annotations/<annotationType>>/<annotationId>/for-<somewhere>/values.bin"
        STORM_LOG_THROW(p.size() == 5 && p[0] == "annotations" && p[4] == "values.bin", storm::exceptions::WrongFormatException,
                        "Unexpected file path '" << srcPath << "'. Expected 'annotations/<annotationType>/<annotationId>/for-<somewhere>/values.bin'.");
        auto const& annotationType = p[1];
        auto const& annotationId = p[2];
        if (annotationType == "aps") {
            auto const& aps = index.annotations.aps;
            STORM_LOG_THROW(aps->contains(annotationId), storm::exceptions::WrongFormatException,
                            "Annotation id '" << annotationId << "' for aps referenced in files but not found in index.");
            auto const annotationType = aps->at(annotationId).type;
            importGenericVector(src, annotationType.value_or(storm::umb::ModelIndex::Annotations::Annotation::Type::Bool), target);
        } else if (annotationType == "rewards") {
            auto const& rewards = index.annotations.rewards;
            STORM_LOG_THROW(rewards->contains(annotationId), storm::exceptions::WrongFormatException,
                            "Annotation id '" << annotationId << "' for rewards referenced in files but not found in index.");
            auto const rewType = rewards->at(annotationId).type;
            STORM_LOG_THROW(rewType.has_value(), storm::exceptions::WrongFormatException,
                            "Reward type is not set in the index, but rewards are present in '" << srcPath << "'.");
            importGenericVector(src, rewType.value(), target);
        }
    }
}

template<typename UmbStructure>
    requires HasFileNames<UmbStructure>
bool importVector(storm::io::ArchiveReadEntry& src, storm::umb::ModelIndex const& index, UmbStructure& umbStructure, std::filesystem::path const& context) {
    static_assert(UmbStructure::FileNames.size() == boost::pfr::tuple_size_v<UmbStructure>, "Number of file names does not match number of fields in struct.");

    // helper to check if src is in the given path
    auto containsSrc = [&src](auto const& other) {
        auto srcPath = getFilePath(src);
        auto srcIt = srcPath.begin();
        for (auto const& o : other) {
            if (srcIt == srcPath.end() || (*srcIt != o && !o.empty())) {
                return false;
            }
            ++srcIt;
        }
        return true;
    };

    // check if the src is in the given context
    if (!containsSrc(context)) {
        return false;
    }
    bool found = false;
    boost::pfr::for_each_field(umbStructure, [&](auto& field, std::size_t i) {
        auto fieldPath = context / std::data(UmbStructure::FileNames)[i];
        // handle the case that we already found the subfield or that this is not the right subfield
        if (found || !containsSrc(fieldPath)) {
            return;
        }

        // load the file into this field, either with a recursive call or directly if the field points to a file
        using FieldType = std::remove_cvref_t<decltype(field)>;
        if constexpr (HasFileNames<FieldType>) {
            found = importVector(src, index, field, fieldPath);
        } else if constexpr (IsOptionalWithFileNames<FieldType>) {
            field.emplace();
            found = importVector(src, index, field.value(), fieldPath);
        } else if constexpr (FileNameMap<FieldType>) {
            for (auto& [key, value] : field) {
                found = importVector(src, index, value, fieldPath / key);
                if (found) {
                    break;
                }
            }
        } else {
            // reaching this point means that we have found the right field
            STORM_LOG_THROW(fieldPath == getFilePath(src), storm::exceptions::WrongFormatException,
                            "Unexpected file paths: " << fieldPath << " != " << getFilePath(src));
            found = true;
            if constexpr (std::is_same_v<FieldType, GenericVector>) {
                importGenericVector(src, index, field);
            } else if constexpr (!std::is_same_v<FieldType, storm::umb::ModelIndex>) {
                importVector(src, field);
            }
        }
    });
    return found;
}

storm::umb::UmbModel fromArchive(std::filesystem::path const& umbArchive, ImportOptions const& /* options */) {
    storm::utility::Stopwatch stopwatch;
    stopwatch.start();
    storm::umb::UmbModel umbModel;
    // First pass: find the index file
    bool indexFound = false;
    for (auto entry : storm::io::openArchive(umbArchive)) {
        if (entry.name() == "index.json") {
            parseIndexFromString(entry.toString(), umbModel.index);
            indexFound = true;
            break;
        }
    }
    STORM_LOG_THROW(indexFound, storm::exceptions::FileIoException, "File 'index.json' not found in UMB archive.");
    STORM_LOG_TRACE("Index file found in umb archive " << umbArchive << ": \n" << storm::dumpJson(storm::json<storm::RationalNumber>(result->index)));
    std::cout << "First pass: Index file loaded in " << stopwatch << " seconds.\n";
    stopwatch.restart();
    // Second pass: load the bin files
    prepareAnnotations(umbModel);
    for (auto entry : storm::io::openArchive(umbArchive)) {
        if (entry.name() == "index.json" || entry.isDir()) {
            continue;  // skip the index file and directories
        }
        bool found = importVector(entry, umbModel.index, umbModel, "");
        STORM_LOG_TRACE("File " << getFilePath(entry) << " in UMB archive " << umbArchive << "' loaded.");
        STORM_LOG_WARN_COND(
            found, "File " << getFilePath(entry) << " in UMB archive " << umbArchive << " will be ignored as it could not be associated with any UMB field.");
    }
    std::cout << "Second pass: bin files loaded in " << stopwatch << " seconds.\n";
    return umbModel;
}

}  // namespace internal

storm::umb::UmbModel importUmb(std::filesystem::path const& umbLocation, ImportOptions const& options) {
    STORM_LOG_THROW(std::filesystem::exists(umbLocation), storm::exceptions::FileIoException, "The given path '" << umbLocation << "' does not exist.");
    return internal::fromArchive(umbLocation, options);
}

}  // namespace storm::umb
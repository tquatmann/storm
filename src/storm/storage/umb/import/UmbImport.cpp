#include "storm/storage/umb/import/UmbImport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/adapters/JsonAdapter.h"
#include "storm/io/ArchiveReader.h"
#include "storm/io/file.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/FileIoException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/WrongFormatException.h"

namespace storm::umb {

namespace internal {

void parseIndexFromDisk(std::filesystem::path const& indexFilePath, storm::umb::ModelIndex& index) {
    storm::json<storm::RationalNumber> parsedStructure;
    std::ifstream file;
    storm::io::openFile(indexFilePath, file);
    file >> parsedStructure;
    storm::io::closeFile(file);
    parsedStructure.get_to(index);
}

void parseIndexFromString(std::string const& indexFileString, storm::umb::ModelIndex& index) {
    std::cout << indexFileString << std::endl;
    storm::json<storm::RationalNumber>::parse(indexFileString).get_to(index);
}

struct UmbPath {
    std::filesystem::path base, filepath;
};

template<typename T>
concept SourceType = std::same_as<T, UmbPath> || std::same_as<T, storm::io::ArchiveReadEntry>;

bool matches(UmbPath const& src, std::filesystem::path const& file) {
    return src.filepath == file;
}

bool matches(storm::io::ArchiveReadEntry const& src, std::filesystem::path const& file) {
    return src.name() == file;
}

std::filesystem::path getFilePath(SourceType auto const& src) {
    if constexpr (std::same_as<decltype(src), UmbPath>) {
        return src.filepath;
    } else {
        return src.name();
    }
}

template<typename T>
bool loadIfMatches(UmbPath const& src, std::filesystem::path const& file, T& objRef) {
    if (matches(src, file)) {
        objRef.emplace(src.base / src.filepath);
        return true;
    }
    return false;
}

template<typename T>
bool loadIfMatches(storm::io::ArchiveReadEntry& src, auto const& file, std::optional<T>& objRef) {
    if (matches(src, file)) {
        if constexpr (std::is_same_v<T, storm::storage::BitVector>) {
            objRef = src.template toVector<bool>();
        } else {
            objRef = src.template toVector<typename T::value_type>();
        }
        return true;
    }
    return false;
}

template<typename EnumType, StorageType Storage>
    requires std::is_enum_v<EnumType>
bool loadGenericVectorIfMatches(SourceType auto& src, std::filesystem::path const& file, EnumType const type, GenericVector<Storage>& objRef) {
    if (!matches(src, file)) {
        return false;
    }
    auto set = [&objRef, &src]<typename T>() {
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(src)>, storm::io::ArchiveReadEntry>) {
            objRef.template set<T>(src.template toVector<T>());
        } else {
            objRef.set(typename GenericVector<Storage>::template Vec<T>(src.base / src.filepath));
        }
    };

    if (type == EnumType::Double) {
        set.template operator()<double>();
    }
    if constexpr (std::is_same_v<EnumType, storm::umb::ModelIndex::Annotation::Type>) {
        if (type == EnumType::Bool) {
            set.template operator()<bool>();
        } else if (type == EnumType::Int32) {
            set.template operator()<int32_t>();
        }
    }
    STORM_LOG_ASSERT(objRef.hasValue(), "Unexpected type with index: " << static_cast<std::underlying_type_t<EnumType>>(type) << ".");
    return true;
}

template<typename EnumType, StorageType Storage>
    requires std::is_enum_v<typename EnumType::E>
bool loadGenericVectorIfMatches(SourceType auto& src, std::filesystem::path const& file, EnumType const type, GenericVector<Storage>& objRef) {
    return loadGenericVectorIfMatches<typename EnumType::E>(src, file, type, objRef);  // convert to inner enum type
}

template<StorageType Storage>
bool loadStatesChoicesBranches(SourceType auto& src, UmbModel<Storage>& umbModel) {
    bool processed{false};

    auto const& ts = umbModel.index.transitionSystem;
    auto& states = umbModel.states;
    processed |= loadIfMatches(src, "state-to-choice.bin", states.stateToChoice);
    processed |= loadIfMatches(src, "state-to-player.bin", states.stateToPlayer);
    processed |= loadIfMatches(src, "initial-states.bin", states.initialStates);
    auto& choices = umbModel.choices;
    processed |= loadIfMatches(src, "choice-to-branch.bin", choices.choiceToBranch);
    processed |= loadIfMatches(src, "choice-to-action.bin", choices.choiceToAction);
    auto& branches = umbModel.branches;
    processed |= loadIfMatches(src, "branch-to-target.bin", branches.branchToTarget);

    using BranchValues = storm::umb::ModelIndex::TransitionSystem::BranchValues;
    using BranchValueType = storm::umb::ModelIndex::TransitionSystem::BranchValueType;
    if (ts.branchValues == BranchValues::Number) {
        STORM_LOG_THROW(ts.branchValueType.has_value(), storm::exceptions::WrongFormatException, "Branch values are numbers, but no type is specified.");
        if (auto const branchVT = *ts.branchValueType; branchVT == BranchValueType::Double) {
            processed |= loadGenericVectorIfMatches(src, "branch-values.bin", branchVT, branches.branchValues);
        } else {
            //  constructFromPathIfExists(umbDir, "branch-to-value.bin", "branch-rational.bin", branches.branchToValue);
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value type " << branchVT << " unhandled.");
        }
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value kind " << ts.branchValues << " unhandled.");
    }
    return processed;
}

template<StorageType Storage>
bool loadAnnotations(SourceType auto& src, UmbModel<Storage>& umbModel) {
    // extract the annotations ID from src
    auto file = getFilePath(src);
    auto fileIt = file.begin();
    auto getNextElementOfPath = [&fileIt, &file]() { return fileIt != file.end() ? *fileIt++ : std::filesystem::path{}; };
    STORM_LOG_ASSERT(getNextElementOfPath() == "annotations", "Unexpected path: " << file);
    std::string const annotationId = getNextElementOfPath();
    auto annotationIndex = umbModel.index.annotations.find(annotationId);
    if (annotationIndex == umbModel.index.annotations.end()) {
        STORM_LOG_WARN("Annotation ID '" << annotationId << "' not found in index file but referenced in file '" << file << "', which will be ignored.");
        return false;
    }
    auto& annotation = umbModel.annotations[annotationId];
    return loadGenericVectorIfMatches(src, "annotations/" + annotationId + "/values.bin", annotationIndex->second.type, annotation.values);
}

std::unique_ptr<UmbModelBase> fromDirectory(std::filesystem::path const& umbDir, ImportOptions const& options) {
    STORM_LOG_THROW(std::filesystem::is_directory(umbDir), storm::exceptions::FileIoException, "The given path is not a directory.");
    auto result = std::make_unique<UmbModel<StorageType::Disk>>();
    internal::parseIndexFromDisk(umbDir / "index.json", result->index);
    std::cout << "umb dir is " << umbDir << "\n";
    for (auto f : std::filesystem::recursive_directory_iterator(umbDir)) {
        std::cout << "Reading file: " << f.path() << " aka " << f << "\n";
        UmbPath entry{umbDir, f.path()};
        if (entry.filepath == "index.json" || f.is_directory() || entry.filepath.empty()) {
            continue;  // skip the index file and directories
                       //        } else if (*entry.filepath.begin() == "annotations") {
                       //            loadAnnotations(entry, *result);
                       //        } else {
                       //            loadStatesChoicesBranches(entry, *result);
        }
    }
    return result;
}

std::unique_ptr<UmbModelBase> fromArchive(std::filesystem::path const& umbArchive, ImportOptions const& options) {
    storm::utility::Stopwatch stopwatch;
    stopwatch.start();
    auto result = std::make_unique<UmbModel<StorageType::Memory>>();
    // First pass: find the index file
    bool indexFound = false;
    for (auto entry : storm::io::openArchive(umbArchive)) {
        if (entry.name() == "index.json") {
            parseIndexFromString(entry.toString(), result->index);
            indexFound = true;
            break;
        }
    }
    STORM_LOG_THROW(indexFound, storm::exceptions::FileIoException, "File 'index.json' not found in UMB archive.");
    std::cout << "First pass: Index file loaded in " << stopwatch << " seconds.\n";
    stopwatch.restart();
    // Second pass: load the bin files
    for (auto entry : storm::io::openArchive(umbArchive)) {
        if (entry.name() == "index.json" || entry.isDir()) {
            continue;  // skip the index file and directories
        } else if (entry.name().starts_with("annotations")) {
            if (!loadAnnotations(entry, *result)) {
                STORM_LOG_WARN("Unable to process file " << entry.name() << " for annotations.");
            }
        } else {
            std::cout << "Reading file: " << entry.name() << "\n";
            if (!loadStatesChoicesBranches(entry, *result)) {
                STORM_LOG_WARN("Unable to process file " << entry.name() << ".");
            }
        }
    }
    std::cout << "Second pass: bin files loaded in " << stopwatch << " seconds.\n";
    return result;
}

}  // namespace internal

std::unique_ptr<UmbModelBase> fromDisk(std::filesystem::path const& umbLocation, ImportOptions const& options) {
    STORM_LOG_THROW(std::filesystem::exists(umbLocation), storm::exceptions::FileIoException, "The given path '" << umbLocation << "' does not exist.");
    if (std::filesystem::is_directory(umbLocation)) {
        return internal::fromDirectory(umbLocation, options);
    } else {
        return internal::fromArchive(umbLocation, options);
    }
}

}  // namespace storm::umb
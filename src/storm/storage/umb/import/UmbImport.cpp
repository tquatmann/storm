#include "storm/storage/umb/import/UmbImport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/io/file.h"
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

template<typename T>
bool constructFromPathIfExists(std::filesystem::path const& file, T& objRef) {
    if (std::filesystem::exists(file)) {
        objRef.emplace(file);
        return true;
    }
    return false;
}

template<typename EnumType>
    requires std::is_enum_v<EnumType>
bool constructGenericVectorFromPathIfExists(std::filesystem::path const& file, EnumType const type, GenericVector<StorageType::Disk>& objRef) {
    if (!std::filesystem::exists(file)) {
        return false;
    }
    if (type == EnumType::Double) {
        objRef.set<double>(typename GenericVector<StorageType::Disk>::Vec<double>(file));
    }
    if constexpr (std::is_same_v<EnumType, storm::umb::ModelIndex::Annotation::Type>) {
        if (type == EnumType::Bool) {
            objRef.set<bool>(typename GenericVector<StorageType::Disk>::Vec<bool>(file));
        } else if (type == EnumType::Int32) {
            objRef.set<int32_t>(typename GenericVector<StorageType::Disk>::Vec<int32_t>(file));
        }
    }
    STORM_LOG_ASSERT(objRef.hasValue(), "Unexpected type with index: " << static_cast<std::underlying_type_t<EnumType>>(type) << ".");
    return true;
}

template<typename EnumType>
    requires std::is_enum_v<typename EnumType::E>
bool constructGenericVectorFromPathIfExists(std::filesystem::path const& file, EnumType const type, GenericVector<StorageType::Disk>& objRef) {
    return constructGenericVectorFromPathIfExists<typename EnumType::E>(file, type, objRef);  // convert to inner enum type
}

void loadStatesChoicesBranchesFromDisk(std::filesystem::path const& umbDir, UmbModel<StorageType::Disk>& umbModel) {
    auto const& ts = umbModel.index.transitionSystem;
    auto& states = umbModel.states;
    constructFromPathIfExists(umbDir / "state-to-choice.bin", states.stateToChoice);
    constructFromPathIfExists(umbDir / "state-to-player.bin", states.stateToPlayer);
    constructFromPathIfExists(umbDir / "initial-states.bin", states.initialStates);
    auto& choices = umbModel.choices;
    constructFromPathIfExists(umbDir / "choice-to-branch.bin", choices.choiceToBranch);
    constructFromPathIfExists(umbDir / "choice-to-action.bin", choices.choiceToAction);
    // constructFromPathIfExists(umbDir, "action-to-action-string.bin", "action-strings.bin", choices.actionStrings);
    auto& branches = umbModel.branches;
    constructFromPathIfExists(umbDir / "branch-to-target.bin", branches.branchToTarget);

    using BranchValues = storm::umb::ModelIndex::TransitionSystem::BranchValues;
    using BranchValueType = storm::umb::ModelIndex::TransitionSystem::BranchValueType;
    if (ts.branchValues == BranchValues::Number) {
        STORM_LOG_THROW(ts.branchValueType.has_value(), storm::exceptions::WrongFormatException, "Branch values are numbers, but no type is specified.");
        if (auto const branchVT = *ts.branchValueType; branchVT == BranchValueType::Double) {
            constructGenericVectorFromPathIfExists(umbDir / "branch-values.bin", branchVT, branches.branchValues);
        } else {
            //  constructFromPathIfExists(umbDir, "branch-to-value.bin", "branch-rational.bin", branches.branchToValue);
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value type " << branchVT << " unhandled.");
        }
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Branch value kind " << ts.branchValues << " unhandled.");
    }
}

void loadAnnotationsFromDisk(std::filesystem::path const& annotationsDir, UmbModel<StorageType::Disk>& umbModel) {
    for (auto const& [annotationId, annotation] : umbModel.index.annotations) {
        std::filesystem::path annotationIdPath{annotationId};
        STORM_LOG_WARN_COND(annotationIdPath.filename() == annotationIdPath, "Unexpeced annotation id " << annotationId);
        auto const dirOfAnnotation = annotationsDir / annotationId;
        STORM_LOG_THROW(std::filesystem::exists(dirOfAnnotation), storm::exceptions::FileIoException,
                        "No files for annotation '" << annotationId << "' referenced in index file.");
        auto& annotationFile = umbModel.annotations[annotationId];
        constructGenericVectorFromPathIfExists(dirOfAnnotation / "values.bin", annotation.type, annotationFile.values);
    }
}

}  // namespace internal

std::unique_ptr<UmbModelBase> fromDisk(std::filesystem::path const& umbDir, ImportOptions const& options) {
    STORM_LOG_THROW(is_directory(umbDir), storm::exceptions::FileIoException, "The given path is not a directory.");
    auto result = std::make_unique<UmbModel<StorageType::Disk>>();
    internal::parseIndexFromDisk(umbDir / "index.json", result->index);
    internal::loadStatesChoicesBranchesFromDisk(umbDir, *result);
    internal::loadAnnotationsFromDisk(umbDir / "annotations", *result);
    return result;
}

}  // namespace storm::umb
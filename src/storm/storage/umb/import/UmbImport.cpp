#include "storm/storage/umb/import/UmbImport.h"

#include "storm/storage/umb/model/UmbModel.h"

#include "storm/io/file.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/FileIoException.h"

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
bool constructFromPathIfExists(std::filesystem::path const& umbDir, std::string_view filename, T& objRef) {
    if (std::filesystem::exists(umbDir / filename)) {
        objRef.emplace(umbDir / filename);
        return true;
    }
    return false;
}

bool constructDoubleVectorFromPathIfExists(std::filesystem::path const& umbDir, std::string_view filename, GenericValueVectorType<StorageType::Disk>& objRef) {
    if (std::filesystem::exists(umbDir / filename)) {
        objRef.set<double>(typename GenericValueVectorType<StorageType::Disk>::Vec<double>(umbDir / filename));
        return true;
    }
    return false;
}

void loadStatesChoicesBranchesFromDisk(std::filesystem::path const& umbDir, UmbModel<StorageType::Disk>& umbModel) {
    auto& states = umbModel.states;
    constructFromPathIfExists(umbDir, "state-to-choice.bin", states.stateToChoice);
    constructFromPathIfExists(umbDir, "state-to-player.bin", states.stateToPlayer);
    constructFromPathIfExists(umbDir, "initial-states.bin", states.initialStates);
    auto& choices = umbModel.choices;
    constructFromPathIfExists(umbDir, "choice-to-branch.bin", choices.choiceToBranch);
    constructFromPathIfExists(umbDir, "choice-to-action.bin", choices.choiceToAction);
    // constructFromPathIfExists(umbDir, "action-to-action-string.bin", "action-strings.bin", choices.actionStrings);
    auto& branches = umbModel.branches;
    constructFromPathIfExists(umbDir, "branch-to-target.bin", branches.branchToTarget);
    // if (umbModel.index.)
    //  constructFromPathIfExists(umbDir, "branch-to-value.bin", "branch-rational.bin", branches.branchToValue);
    constructDoubleVectorFromPathIfExists(umbDir, "branch-to-value.bin", branches.branchToValue);
}

}  // namespace internal

std::unique_ptr<UmbModelBase> fromDisk(std::filesystem::path const& umbDir, ImportOptions const& options) {
    STORM_LOG_THROW(is_directory(umbDir), storm::exceptions::FileIoException, "The given path is not a directory.");
    auto result = std::make_unique<UmbModel<StorageType::Disk>>();
    internal::parseIndexFromDisk(umbDir / "index.json", result->index);
    internal::loadStatesChoicesBranchesFromDisk(umbDir, *result);
    return result;
}

}  // namespace storm::umb
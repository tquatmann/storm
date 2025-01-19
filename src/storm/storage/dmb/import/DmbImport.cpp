#include "storm/storage/dmb/import/DmbImport.h"

#include "storm/storage/dmb/model/DmbModel.h"

#include "storm/io/file.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/FileIoException.h"

namespace storm::dmb {

namespace internal {

void parseIndexFromDisk(std::filesystem::path const& indexFilePath, storm::dmb::ModelIndex& index) {
    storm::json<storm::RationalNumber> parsedStructure;
    std::ifstream file;
    storm::io::openFile(indexFilePath, file);
    file >> parsedStructure;
    storm::io::closeFile(file);
    parsedStructure.get_to(index);
}

template<typename T>
bool constructFromPathIfExists(std::filesystem::path const& dmbDir, std::string_view filename, T& objRef) {
    if (std::filesystem::exists(dmbDir / filename)) {
        objRef.emplace(dmbDir / filename);
        return true;
    }
    return false;
}

bool constructDoubleVectorFromPathIfExists(std::filesystem::path const& dmbDir, std::string_view filename, GenericValueVectorType<StorageType::Disk>& objRef) {
    if (std::filesystem::exists(dmbDir / filename)) {
        objRef.set<double>(typename GenericValueVectorType<StorageType::Disk>::Vec<double>(dmbDir / filename));
        return true;
    }
    return false;
}

void loadStatesChoicesBranchesFromDisk(std::filesystem::path const& dmbDir, DmbModel<StorageType::Disk>& dmbModel) {
    auto& states = dmbModel.states;
    constructFromPathIfExists(dmbDir, "state-to-choice.bin", states.stateToChoice);
    constructFromPathIfExists(dmbDir, "state-to-player.bin", states.stateToPlayer);
    // constructFromPathIfExists(dmbDir, "initial-states.bin", states.initialStates);
    auto& choices = dmbModel.choices;
    constructFromPathIfExists(dmbDir, "choice-to-branch.bin", choices.choiceToBranch);
    constructFromPathIfExists(dmbDir, "choice-to-action.bin", choices.choiceToAction);
    // constructFromPathIfExists(dmbDir, "action-to-action-string.bin", "action-strings.bin", choices.actionStrings);
    auto& branches = dmbModel.branches;
    constructFromPathIfExists(dmbDir, "branch-to-target.bin", branches.branchToTarget);
    // if (dmbModel.index.)
    //  constructFromPathIfExists(dmbDir, "branch-to-value.bin", "branch-rational.bin", branches.branchToValue);
    constructDoubleVectorFromPathIfExists(dmbDir, "branch-to-value.bin", branches.branchToValue);
}

}  // namespace internal

std::unique_ptr<DmbModelBase> fromDisk(std::filesystem::path const& dmbDir, ImportOptions const& options) {
    STORM_LOG_THROW(is_directory(dmbDir), storm::exceptions::FileIoException, "The given path is not a directory.");
    auto result = std::make_unique<DmbModel<StorageType::Disk>>();
    internal::parseIndexFromDisk(dmbDir / "index.json", result->index);
    internal::loadStatesChoicesBranchesFromDisk(dmbDir, *result);
    return result;
}

}  // namespace storm::dmb
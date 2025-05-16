#include "storm/storage/umb/model/UmbModel.h"

namespace storm::umb {

bool validateCsr(auto&& csr, uint64_t numMappedElements, uint64_t expectedlastEntry, bool silent) {
    if (csr) {
        // Check if csr has expected length and the form {0, ..., expectedlastEntry}
        if (csr->size() != numMappedElements + 1) {
            STORM_LOG_ERROR_COND(silent, "CSR has unexpected size: " << csr->size() << " != " << (numMappedElements + 1));
            return false;
        }
        if (csr.value()[0] != 0) {
            STORM_LOG_ERROR_COND(silent, "CSR has unexpected first entry: " << csr.value()[0] << " != 0");
            return false;
        }
        if (csr.value()[numMappedElements] != expectedlastEntry) {
            STORM_LOG_ERROR_COND(silent, "CSR has unexpected last entry: " << csr.value()[numMappedElements] << " != " << expectedlastEntry);
            return false;
        }
    } else if (numMappedElements != expectedlastEntry) {  // we assume a 1:1 mapping
        STORM_LOG_ERROR_COND(silent, "CSR is not given and the default 1:1 mapping {0, ... ,"
                                         << numMappedElements << "} does not match. Expected the mapping to end with '" << expectedlastEntry << "'");
        return false;
    }
    return true;
}

template<StorageType Storage>
bool UmbModel<Storage>::validate(bool silent) const {
    auto const& tsIndex = index.transitionSystem;
    bool isValid = true;

    // validate counts
    auto checkNum = [&silent](uint64_t num, auto&& name) {
        if (num == storm::umb::ModelIndex::TransitionSystem::InvalidNumber) {
            STORM_LOG_ERROR_COND(silent, "Number of " << name << " is not set.");
            return false;
        }
        return true;
    };
    isValid &= checkNum(tsIndex.numPlayers, "players");
    isValid &= checkNum(tsIndex.numStates, "states");
    isValid &= checkNum(tsIndex.numInitialStates, "initial-states");
    isValid &= checkNum(tsIndex.numChoices, "choices");
    isValid &= checkNum(tsIndex.numActions, "actions");
    isValid &= checkNum(tsIndex.numBranches, "branches");

    // Validate CSR mappings
    auto sizeOr = [](auto&& input, std::size_t defaultValue) { return input.has_value() ? input->size() : defaultValue; };
    if (!validateCsr(states.stateToChoice, tsIndex.numStates, tsIndex.numChoices, silent)) {
        STORM_LOG_ERROR_COND(silent, "State to choice mapping is invalid.");
        isValid = false;
    }
    if (!validateCsr(choices.choiceToBranch, tsIndex.numChoices, tsIndex.numBranches, silent)) {
        STORM_LOG_ERROR_COND(silent, "Choice to branch mapping is invalid.");
        isValid = false;
    }
    // todo: check what the default behavior of action strings is.
    //    if (!validateCsr(choices.actionToActionString, tsIndex.numActions, sizeOr(choices.actionStrings, 0), silent)) {
    //        STORM_LOG_ERROR_COND(silent, "Action string mapping is invalid.");
    //        isValid = false;
    //    }

    // Validate other inputs
    if (!branches.branchToTarget.has_value() || branches.branchToTarget->size() != tsIndex.numBranches) {
        STORM_LOG_ERROR_COND(silent, "Branch to target mapping is missing or has invalid size: " << sizeOr(branches.branchToTarget, 0));
        isValid = false;
    }
    if (!branches.branchProbabilities.hasValue()) {
        STORM_LOG_ERROR_COND(silent, "Branch probabilities are missing.");
        isValid = false;
    }

    // TODO: add more validations
    // If prob type is rational, probabilities are either rational or uint64
    // if annotation forStates/forChoices/forBranches is set, values are set too
    // annotations in index are unique, and present as files
    return isValid;
}

template class UmbModel<StorageType::Disk>;
template class UmbModel<StorageType::Memory>;

}  // namespace storm::umb

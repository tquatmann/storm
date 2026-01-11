#include "storm/storage/umb/model/Validation.h"

#include <sstream>
#include <string_view>

#include "storm/storage/umb/model/UmbModel.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/WrongFormatException.h"

namespace storm::umb {

namespace validation {
bool validateCsr(auto const& csr, std::string_view const name, uint64_t numMappedElements, uint64_t expectedlastEntry, std::ostream& err) {
    std::stringstream err_reason;
    if (csr) {
        // Check if csr has expected length and the form {0, ..., expectedlastEntry}
        if (csr->size() != numMappedElements + 1) {
            err_reason << "CSR has unexpected size: " << csr->size() << " != " << (numMappedElements + 1) << ".";
        }
        if (csr.value()[0] != 0) {
            err_reason << "CSR has unexpected first entry: " << csr.value()[0] << " != 0" << ".";
        }
        if (csr.value()[numMappedElements] != expectedlastEntry) {
            err_reason << "CSR has unexpected last entry: " << csr.value()[numMappedElements] << " != " << expectedlastEntry << ".";
        }
    } else if (numMappedElements != expectedlastEntry) {  // we assume a 1:1 mapping
        err_reason << "CSR is not given and the default 1:1 mapping {0, ... ," << numMappedElements << "} does not match. Expected the mapping to end with '"
                   << expectedlastEntry << "'" << ".";
    }
    if (!err_reason.view().empty()) {
        err << "Validation error in CSR mapping '" << name << "':\n\t" << err_reason.str() << "\n";
        return false;
    }
    return true;
}
}  // namespace validation

bool validate(storm::umb::UmbModel const& umbModel, std::ostream& err) {
    auto const& index = umbModel.index;
    auto const& tsIndex = index.transitionSystem;
    bool isValid = true;

    // validate counts
    auto checkNum = [&err](uint64_t num, auto&& name) {
        if (num == storm::umb::ModelIndex::TransitionSystem::InvalidNumber) {
            err << "Number of " << name << " is not set.\n";
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
    isValid &= validation::validateCsr(umbModel.states.stateToChoice, "state-to-choice", tsIndex.numStates, tsIndex.numChoices, err);
    isValid &= validation::validateCsr(umbModel.choices.choiceToBranch, "choice-to-branch", tsIndex.numChoices, tsIndex.numBranches, err);
    // todo: check what the default behavior of action strings is.
    // isValid &= validation::validateCsr(umbModel.choices.actionToActionString, "action-to-action-strings", tsIndex.numActions,
    //                                   sizeOr(umbModel.choices.actionStrings, 0), err);

    // Validate state valuations
    if (index.stateValuations.has_value() != umbModel.states.stateValuations.has_value()) {
        if (index.stateValuations.has_value()) {
            err << "State valuations described in index file but not present.\n";
        } else {
            err << "State valuations present but not described in index file.\n";
        }
        isValid = false;
    } else if (index.stateValuations.has_value()) {
        uint64_t const alignment = index.stateValuations->alignment;
        uint64_t const numBytes = umbModel.states.stateValuations->size();
        if (alignment == 0) {
            err << "State valuation alignment is 0.\n";
            isValid = false;
        }
        if (numBytes % alignment != 0) {
            err << "State valuation data size is " << umbModel.states.stateValuations->size() << " which is not a multiple of the alignment '" << alignment
                << "'.\n";
            isValid = false;
        }
        isValid &= validation::validateCsr(umbModel.states.stateToValuation, "state-to-valuation", tsIndex.numStates, numBytes / alignment, err);
    }

    // Validate other inputs
    if (!umbModel.branches.branchToTarget.has_value()) {
        err << "Branch to target mapping is missing.\n";
        isValid = false;
    } else if (umbModel.branches.branchToTarget->size() != tsIndex.numBranches) {
        err << "Branch to target mapping has invalid size: " << umbModel.branches.branchToTarget->size() << " != " << tsIndex.numBranches << ".\n";
        isValid = false;
    }
    if (!umbModel.branches.branchProbabilities.hasValue()) {
        err << "Branch probabilities are missing.";
        isValid = false;
    }

    // TODO: add more validations
    // If prob type is rational, probabilities are either rational or uint64
    // if annotation forStates/forChoices/forBranches is set, values are set too
    // annotations in index are unique, and present as files
    // rewards are numeric, aps are boolean
    // annotations appliesTo is consistent with present files, annotation types are consistent (aps = bool, rewards = numeric)
    return isValid;
}

void validateOrThrow(storm::umb::UmbModel const& umbModel) {
    std::stringstream errors;
    if (!validate(umbModel, errors)) {
        STORM_LOG_THROW(false, storm::exceptions::WrongFormatException,
                        "UMB model " << umbModel.getShortModelInformation() << " is invalid:\n"
                                     << errors.str());
    }
}

}  // namespace storm::umb

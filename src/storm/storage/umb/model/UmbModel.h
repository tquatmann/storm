#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "storm/storage/BitVector.h"
#include "storm/storage/umb/model/GenericVector.h"
#include "storm/storage/umb/model/ModelIndex.h"
#include "storm/storage/umb/model/Types.h"

namespace storm::umb {
class UmbModel {
   public:
    ModelIndex index;

    struct States {
        CSR stateToChoice;
        TO1<uint32_t> stateToPlayer;
        TO1<bool> initialStates;
        TO1<bool> markovianStates;
        CSR stateToExitRate;
        SEQ<AnyValueType> exitRates;
        CSR stateToValuation;
        SEQ<char> stateValuations;
        auto static constexpr FileNames = {"state-to-choice.bin",    "state-to-player.bin", "initial-states.bin",     "markovian-states.bin",
                                           "state-to-exit-rate.bin", "exit-rates.bin",      "state-to-valuation.bin", "state-valuations.bin"};
    } states;
    struct Choices {
        CSR choiceToBranch;
        TO1<uint32_t> choiceToAction;
        CSR actionToActionString;
        SEQ<char> actionStrings;
        auto static constexpr FileNames = {"choice-to-branch.bin", "choice-to-action.bin", "action-to-action-string.bin", "action-strings.bin"};
    } choices;
    struct Branches {
        TO1<uint64_t> branchToTarget;
        CSR branchToProbability;
        SEQ<AnyValueType> branchProbabilities;
        auto static constexpr FileNames = {"branch-to-target.bin", "branch-to-probability.bin", "branch-probabilities.bin"};
    } branches;

    struct Annotation {
        struct Values {
            CSR toValue;
            SEQ<AnyValueType> values;
            auto static constexpr FileNames = {"to-value.bin", "values.bin"};
        };
        std::optional<Values> forStates, forChoices, forBranches;
        auto static constexpr FileNames = {"for-states/", "for-choices/", "for-branches/"};
    };
    std::map<std::string, Annotation> rewards, aps;

    auto static constexpr FileNames = {"index.json", "", "", "", "annotations/rewards/", "annotations/aps/"};

    /*!
     * Retrieves a short string that can be used to refer to the model in user output.
     */
    std::string getShortModelInformation() const;

    /*!
     * Retrieves a string that describes the model in more detail.
     */
    std::string getModelInformation() const;

    /*!
     * Validates the given UMB model and writes potential errors to the given output stream.
     * @return true if the UMB model is valid.
     */
    bool validate(std::ostream& errors) const;

    /*!
     * Validates the UmbModel. If it is invalid, an exception is thrown.
     */
    void validateOrThrow() const;
};

}  // namespace storm::umb
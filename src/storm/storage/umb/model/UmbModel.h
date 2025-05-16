#pragma once

#include <memory>

#include "storm/storage/umb/model/GenericVector.h"
#include "storm/storage/umb/model/ModelIndex.h"
#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/UmbModelBase.h"
#include "storm/storage/umb/model/VectorType.h"

namespace storm::umb {

template<StorageType Storage>
class UmbModel {
   public:
    template<typename T>
    using OptionalVec = std::optional<VectorType<T, Storage>>;
    struct AnyType {};

    template<typename T>
    struct TO1Helper {
        using type = OptionalVec<T>;
    };
    template<>
    struct TO1Helper<AnyType> {
        using type = GenericVector<Storage>;
    };
    template<typename T>
    using TO1 = typename TO1Helper<T>::type;

    template<typename T>
    using SEQ = TO1<T>;
    using CSR = OptionalVec<uint64_t>;

    ModelIndex index;

    struct States {
        CSR stateToChoice;
        TO1<uint32_t> stateToPlayer;
        TO1<bool> initialStates;
        auto static constexpr FileNames = {"state-to-choice.bin", "state-to-player.bin", "initial-states.bin"};
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
        SEQ<AnyType> branchProbabilities;
        auto static constexpr FileNames = {"branch-to-target.bin", "branch-to-probability.bin", "branch-probabilities.bin"};
    } branches;

    struct Annotation {
        struct Values {
            CSR toValue;
            SEQ<AnyType> values;
            auto static constexpr FileNames = {"to-value.bin", "values.bin"};
            // TODO: Support for distributions, string-valued annotations, ...
        };
        std::optional<Values> forStates, forChoices, forBranches;
        auto static constexpr FileNames = {"for-states/", "for-choices/", "for-branches/"};
    };
    std::map<std::string, Annotation> rewards;

    struct BooleanAnnotation {
        struct Values {
            TO1<bool> values;
            auto static constexpr FileNames = {"values.bin"};
        };
        std::optional<Values> forStates, forChoices, forBranches;
        auto static constexpr FileNames = {"for-states/", "for-choices/", "for-branches/"};
    };
    std::map<std::string, BooleanAnnotation> atomicPropositions;

    auto static constexpr FileNames = {"index.json", "", "", "", "annotations/rewards/", "annotations/aps/"};

    bool validate(bool silent = false) const;
};

}  // namespace storm::umb
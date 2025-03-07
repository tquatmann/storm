#pragma once

#include <memory>

#include "storm/storage/umb/model/GenericValueVector.h"
#include "storm/storage/umb/model/ModelIndex.h"
#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/UmbModelBase.h"
#include "storm/storage/umb/model/VectorType.h"

namespace storm::umb {

template<StorageType Storage>
class UmbModel : public UmbModelBase {
   public:
    template<typename T>
    using OptionalVec = std::optional<VectorType<T, Storage>>;
    struct GenericValue {};

    template<typename T>
    struct TO1Helper {
        using type = OptionalVec<T>;
    };
    template<>
    struct TO1Helper<GenericValue> {
        using type = GenericValueVectorType<Storage>;
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
        SEQ<bool> initialStates;
    } states;
    struct Choices {
        CSR choiceToBranch;
        TO1<uint32_t> choiceToAction;
        SEQ<std::string> actionStrings;
    } choices;
    struct Branches {
        TO1<uint64_t> branchToTarget;
        TO1<GenericValue> branchToValue;
    } branches;

    virtual bool isStorageType(StorageType storageType) const override;
    virtual ModelIndex const& getIndex() const override;
};

}  // namespace storm::umb
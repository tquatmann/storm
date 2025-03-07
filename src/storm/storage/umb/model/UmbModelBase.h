#pragma once

#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {
class UmbModelBase {
   public:
    virtual ~UmbModelBase() = default;

    virtual bool isStorageType(StorageType storageType) const = 0;
    virtual ModelIndex const& getIndex() const = 0;
    template<StorageType Storage>
    UmbModel<Storage> const& as() const {
        return dynamic_cast<UmbModel<Storage> const&>(*this);
    }
    template<StorageType Storage>
        requires(Storage == StorageType::Memory)
    UmbModel<Storage>& as() {
        return dynamic_cast<UmbModel<Storage>&>(*this);
    }
    // TODO: to_memory
};

}  // namespace storm::umb

#pragma once

#include "storm/storage/dmb/model/DmbModelForward.h"
#include "storm/storage/dmb/model/StorageType.h"

namespace storm::dmb {
class DmbModelBase {
   public:
    virtual ~DmbModelBase() = default;

    virtual bool isStorageType(StorageType storageType) const = 0;
    virtual ModelIndex const& getIndex() const = 0;
    template<StorageType Storage>
    DmbModel<Storage> const& as() const {
        return dynamic_cast<DmbModel<Storage> const&>(*this);
    }
    template<StorageType Storage>
        requires(Storage == StorageType::Memory)
    DmbModel<Storage>& as() {
        return dynamic_cast<DmbModel<Storage>&>(*this);
    }
    // TODO: to_memory
};

}  // namespace storm::dmb

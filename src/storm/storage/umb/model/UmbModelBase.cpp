#include "storm/storage/umb/model/UmbModelBase.h"

#include "storm/storage/umb/model/UmbModel.h"

namespace storm::umb {

bool UmbModelBase::isStorageType(StorageType storageType) const {
    switch (storageType) {
        case StorageType::Disk:
            return std::holds_alternative<UmbModelPtr<StorageType::Disk>>(data);
        case StorageType::Memory:
            return std::holds_alternative<UmbModelPtr<StorageType::Memory>>(data);
        default:
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected storage type.");
    }
}

ModelIndex const& UmbModelBase::getIndex() const {
    if (isStorageType(StorageType::Disk)) {
        return as<StorageType::Disk>().index;
    } else {
        return as<StorageType::Memory>().index;
    }
}

bool UmbModelBase::validate(bool silent) const {
    if (isStorageType(StorageType::Disk)) {
        return as<StorageType::Disk>().validate(silent);
    } else if (isStorageType(StorageType::Memory)) {
        return as<StorageType::Memory>().validate(silent);
    } else {
        return false;
    }
}

template<StorageType Storage>
UmbModel<Storage> const& UmbModelBase::as() const {
    return *std::get<UmbModelPtr<Storage>>(data);
}

template<StorageType Storage>
    requires(Storage == StorageType::Memory)
UmbModel<Storage>& UmbModelBase::as() {
    return *std::get<UmbModelPtr<Storage>>(data);
}

template UmbModel<StorageType::Disk> const& UmbModelBase::as<StorageType::Disk>() const;
template UmbModel<StorageType::Memory> const& UmbModelBase::as<StorageType::Memory>() const;
template UmbModel<StorageType::Memory>& UmbModelBase::as<StorageType::Memory>();

}  // namespace storm::umb

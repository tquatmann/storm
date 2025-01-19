#include "storm/storage/dmb/model/DmbModel.h"

namespace storm::dmb {

template<StorageType Storage>
bool DmbModel<Storage>::isStorageType(StorageType storageType) const {
    return storageType == Storage;
}

template<StorageType Storage>
ModelIndex const& DmbModel<Storage>::getIndex() const {
    return index;
}

template class DmbModel<StorageType::Disk>;
template class DmbModel<StorageType::Memory>;

}  // namespace storm::dmb

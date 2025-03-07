#include "storm/storage/umb/model/UmbModel.h"

namespace storm::umb {

template<StorageType Storage>
bool UmbModel<Storage>::isStorageType(StorageType storageType) const {
    return storageType == Storage;
}

template<StorageType Storage>
ModelIndex const& UmbModel<Storage>::getIndex() const {
    return index;
}

template class UmbModel<StorageType::Disk>;
template class UmbModel<StorageType::Memory>;

}  // namespace storm::umb

#pragma once

#include "storm/adapters/RationalNumberForward.h"
#include "storm/io/BinaryFileViewer.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/UmbBitVector.h"

namespace storm::umb {

namespace internal {
template<typename T, StorageType Storage>
struct VectorTypeHelper;

template<>
struct VectorTypeHelper<storm::RationalNumber, StorageType::Disk> {
    using type = std::vector<storm::RationalNumber>;  // todo
};

template<>
struct VectorTypeHelper<std::string, StorageType::Disk> {
    using type = std::vector<std::string>;  // todo
};

template<>
struct VectorTypeHelper<bool, StorageType::Disk> {
    using type = UmbBitVector<StorageType::Disk>;
};

template<>
struct VectorTypeHelper<bool, StorageType::Memory> {
    using type = storm::storage::BitVector;
};

template<typename T>
struct VectorTypeHelper<T, StorageType::Disk> {
    using type = storm::io::BinaryFileViewer<T, std::endian::little>;
};

template<typename T>
struct VectorTypeHelper<T, StorageType::Memory> {
    using type = std::vector<T>;
};

}  // namespace internal

template<typename T, StorageType Storage>
using VectorType = typename internal::VectorTypeHelper<T, Storage>::type;
}  // namespace storm::umb
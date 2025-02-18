#pragma once
#include <variant>

#include "storm/io/BinaryFileViewer.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/dmb/model/StorageType.h"

namespace storm::dmb {

template<StorageType Storage>
    requires(Storage != StorageType::Memory)
class DmbBitVector {
   public:
    using UnderlyingType = storm::io::BinaryFileViewer<uint64_t, std::endian::little>;

    DmbBitVector(std::filesystem::path const& path)
        requires(Storage == StorageType::Disk)
        : data(path) {
        // Intentionally left empty
    }
    
    storm::storage::BitVector getAsBitVectorAutoSize() const {
        return getAsBitVector(data.size() * 64ull);
    }

    storm::storage::BitVector getAsBitVector(uint64_t size) const {
        STORM_LOG_ASSERT(size <= data.size() * 64ull, "Invalid size. expected size=" << size << " but data has size=" << data.size() * 64ull << ".");
        storm::storage::BitVector result(size, false);
        for (uint64_t bucketIndex = 0; auto bits : data) {
            result.setBucket(bucketIndex, bits);
            ++bucketIndex;
        }
        return result;
    }

   private:
    UnderlyingType data;
};
}  // namespace storm::dmb
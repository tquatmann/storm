#pragma once
#include <variant>

#include "storm/io/BinaryFileViewer.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/umb/model/StorageType.h"
#include "storm/utility/bitoperations.h"

namespace storm::umb {

class UmbBitVector {
   public:
    using UnderlyingType = storm::io::BinaryFileViewer<uint64_t, std::endian::little>;

    UmbBitVector(std::filesystem::path const& path) : data(path) {
        // Intentionally left empty
    }

    storm::storage::BitVector getAsBitVectorAutoSize() const {
        return getAsBitVector(data.size() * 64ull);
    }

    storm::storage::BitVector getAsBitVector(uint64_t size) const {
        STORM_LOG_ASSERT(size <= data.size() * 64ull, "Invalid size. expected size=" << size << " but data has size=" << data.size() * 64ull << ".");
        storm::storage::BitVector result(size, false);
        for (uint64_t bucketIndex = 0; auto bits : data) {
            result.setBucket(bucketIndex,
                             storm::utility::reverseBits(
                                 bits));  // Our bit vectors store the items in reverse order, i.e., the first item is indicated by the most significant bit
            ++bucketIndex;
        }
        return result;
    }

   private:
    UnderlyingType data;
};
}  // namespace storm::umb
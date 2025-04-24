#include "storm/storage/bisimulation/Partition.h"

#include <ranges>

namespace storm::storage::bisimulation {

Partition::Partition(ElementIndex numElements) : blockIndices(numElements, false), elementToBlockIndex(numElements, 0) {
    auto indexRange = std::ranges::iota_view<ElementIndex, ElementIndex>(0, numElements);
    blockContents.assign(indexRange.begin(), indexRange.end());
    blockIndices.set(0);
}

std::size_t Partition::getNumberOfBlocks() const {
    return blockIndices.getNumberOfSetBits();
}

std::size_t Partition::getNumberOfElements() const {
    return elementToBlockIndex.size();
}

typename Partition::Block Partition::getBlockOfElement(ElementIndex element) const {
    return getBlockFromIndex(elementToBlockIndex[element]);
}

bool Partition::contains(Block const& block, ElementIndex element) const {
    STORM_LOG_ASSERT(element < getNumberOfElements(), "Invalid element index '" << element << "'.");
    auto const eBlockIndex = elementToBlockIndex[element];
    auto const blockStart = getBlockIndex(block);
    return eBlockIndex >= blockStart && eBlockIndex < blockStart + block.size();
}

bool Partition::checkBlockValidity(const storm::storage::bisimulation::Partition::Block& block) const {
    if (block.empty()) {
        STORM_LOG_ERROR("Block is empty.");
        return false;
    }
    if (block.data() < blockContents.data() || block.data() + block.size() > blockContents.data() + blockContents.size()) {
        STORM_LOG_ERROR("Block does not belong to this partition. Accidentally copy of partition?");
        return false;
    }
    auto const blockIndex = std::distance(blockContents.data(), block.data());  // can't use ::getBlockIndex() to avoid infinite recursion
    auto const blockEndIndex = blockIndex + block.size();
    if (!blockIndices.get(blockIndex) || (blockEndIndex < blockIndices.size() && !blockIndices.get(blockEndIndex))) {
        STORM_LOG_ERROR("Block indices are not set correctly.");
        return false;
    }
    if (elementToBlockIndex[block[0]] != blockIndex) {
        STORM_LOG_ERROR("Block index of first element is not set correctly.");
        return false;
    }
    if (std::any_of(block.begin(), block.end(), [this, blockIndex, blockEndIndex](ElementIndex const& e) {
            return elementToBlockIndex[e] < blockIndex || elementToBlockIndex[e] > blockEndIndex;
        })) {
        STORM_LOG_ERROR("Block index of one of its elements is not set correctly.");
        return false;
    }
    return true;
}

bool Partition::isProperSuperBlock(Block const& block) const {
    auto const blockStart = getBlockIndex(block);
    auto const blockEnd = blockStart + block.size();
    return blockIndices.getNextSetIndex(blockStart + 1) != blockEnd;
}

typename Partition::BlockIndex Partition::getBlockIndex(Block const& block) const {
    STORM_LOG_ASSERT(checkBlockValidity(block), "Tried to get the index of an invalid block.");
    return std::distance(blockContents.data(), block.data());
}

typename Partition::Block Partition::getBlockFromIndexRange(BlockIndex const start, BlockIndex const end) const {
    STORM_LOG_ASSERT(start < end && end <= blockIndices.size() && blockIndices.get(start), "Invalid block index range [" << start << ", " << end << ").");
    STORM_LOG_ASSERT(end == blockIndices.size() || blockIndices.get(end), "Invalid block index end in range [" << start << ", " << end << ").");
    return {blockContents.begin() + start, blockContents.begin() + end};
}

typename Partition::Block Partition::getBlockFromIndex(BlockIndex const index) const {
    STORM_LOG_ASSERT(index < blockIndices.size() && blockIndices.get(index), "Invalid block index " << index << ".");
    auto const blockEnd = blockIndices.getNextSetIndex(index + 1);
    return getBlockFromIndexRange(index, blockEnd);
}

std::ostream& operator<<(std::ostream& os, const Partition& partition) {
    os << "Partition (" << partition.getNumberOfBlocks() << " block(s), " << partition.getNumberOfElements() << " element(s)): {\n";
    partition.forEachBlock([&os](Partition::Block const& block) {
        os << "\t{";
        for (bool first{true}; auto const e : block) {
            if (first) {
                first = false;
            } else {
                os << ", ";
            }
            os << e;
        }
        os << "}\n";
    });
    os << "}";
    return os;
}

}  // namespace storm::storage::bisimulation
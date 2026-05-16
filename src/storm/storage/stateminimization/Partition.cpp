#include "storm/storage/stateminimization/Partition.h"

namespace storm::storage::stateminimization {

Partition::Partition(ElementIndex numElements)
    : blockIndices(numElements, false), elementToBlockIndex(numElements, 0), blockEndCache(numElements, std::numeric_limits<BlockIndex>::max()) {
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

Partition::Block Partition::getBlockOfElement(ElementIndex element) const {
    return getBlockFromIndex(elementToBlockIndex[element]);
}

bool Partition::contains(Block const& block, ElementIndex element) const {
    STORM_LOG_ASSERT(element < getNumberOfElements(), "Invalid element index '" << element << "'.");
    auto const eBlockIndex = elementToBlockIndex[element];
    auto const blockStart = getBlockIndex(block);
    return eBlockIndex >= blockStart && eBlockIndex < blockStart + block.size();
}

Partition::Block Partition::getBlockFromIndex(BlockIndex index) const {
    STORM_LOG_ASSERT(index < blockIndices.size() && blockIndices.get(index), "Invalid block index " << index << ".");

    if (blockEndCache.size() <= index) {
        blockEndCache.resize(blockIndices.size(), std::numeric_limits<BlockIndex>::max());
    }

    if (blockEndCache[index] == std::numeric_limits<BlockIndex>::max()) {
        blockEndCache[index] = blockIndices.getNextSetIndex(index + 1);
    }

    return getBlockFromIndexRange(index, blockEndCache[index]);
}

Partition::Block Partition::getBlockFromIndexRange(BlockIndex start, BlockIndex end) const {
    STORM_LOG_ASSERT(start < end && end <= blockIndices.size() && blockIndices.get(start), "Invalid block index range [" << start << ", " << end << ").");
    STORM_LOG_ASSERT(end == blockIndices.size() || blockIndices.get(end), "Invalid block index end in range [" << start << ", " << end << ").");
    return {blockContents.begin() + start, blockContents.begin() + end};
}

Partition::BlockIndex Partition::getBlockIndex(Block const& block) const {
    STORM_LOG_ASSERT(checkBlockValidity(block), "Tried to get the index of an invalid block.");
    return std::distance(blockContents.data(), block.data());
}

void Partition::invalidateCache(Partition::BlockIndex index) {
    if (index < blockEndCache.size()) {
        blockEndCache[index] = std::numeric_limits<BlockIndex>::max();
    }
}

bool Partition::checkBlockValidity(Block const& block) const {
    if (block.empty())
        return false;
    if (block.data() < blockContents.data() || block.data() + block.size() > blockContents.data() + blockContents.size())
        return false;

    auto const blockIndex = std::distance(blockContents.data(), block.data());
    auto const blockEndIndex = blockIndex + block.size();

    if (!blockIndices.get(blockIndex) || (blockEndIndex < blockIndices.size() && !blockIndices.get(blockEndIndex)))
        return false;
    if (elementToBlockIndex[block[0]] != blockIndex)
        return false;

    return !std::any_of(block.begin(), block.end(), [this, blockIndex, blockEndIndex](ElementIndex const& e) {
        return elementToBlockIndex[e] < blockIndex || elementToBlockIndex[e] > blockEndIndex;
    });
}

bool Partition::isProperSuperBlock(Block const& block) const {
    auto const blockStart = getBlockIndex(block);
    auto const blockEnd = blockStart + block.size();
    return blockIndices.getNextSetIndex(blockStart + 1) != blockEnd;
}

std::ostream& operator<<(std::ostream& os, const Partition& partition) {
    os << "Partition (" << partition.getNumberOfBlocks() << " block(s), " << partition.getNumberOfElements() << " element(s)): {\n";
    partition.forEachBlock([&os](Partition::Block const& block) {
        os << "\t{";
        for (bool first = true; auto const e : block) {
            if (!first)
                os << ", ";
            first = false;
            os << e;
        }
        os << "}\n";
    });
    os << "}";
    return os;
}

}  // namespace storm::storage::stateminimization

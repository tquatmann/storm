#include "storm/storage/stateminimization/Partition.h"

namespace storm::storage::stateminimization {

Partition::NonSuperBlockSet::NonSuperBlockSet(Partition const& partition) : partition(partition), blockIndices(partition.getNumberOfElements(), false) {}

std::size_t Partition::NonSuperBlockSet::size() const {
    return containedBlockIndices.size();
}

bool Partition::NonSuperBlockSet::empty() const {
    return containedBlockIndices.empty();
}

bool Partition::NonSuperBlockSet::contains(Block const& block) const {
    return blockIndices.get(partition.getBlockIndex(block));
}

void Partition::NonSuperBlockSet::insert(Block const& block) {
    STORM_LOG_ASSERT(!block.empty(), "Cannot insert empty block.");
    STORM_LOG_ASSERT(!partition.isProperSuperBlock(block), "Cannot insert block that is a proper superblock.");
    if (BlockIndex i = partition.getBlockIndex(block); !blockIndices.get(i)) {
        blockIndices.set(i, true);
        containedBlockIndices.push_back(i);
    }
}

Partition::Block Partition::NonSuperBlockSet::pop() {
    STORM_LOG_ASSERT(!empty(), "Cannot pop from empty blockset.");
    auto const i = containedBlockIndices.back();
    blockIndices.set(i, false);
    containedBlockIndices.pop_back();
    return partition.getBlockFromIndex(i);
}

Partition::Partition(ElementIndex numElements) : numBlocks(1) {
    auto indexRange = std::ranges::iota_view<ElementIndex, ElementIndex>(0, numElements);
    blockContents.assign(indexRange.begin(), indexRange.end());
    blockContentsInverse = blockContents;
    elementToBlockIndex.assign(numElements, 0);
    blockSizes.assign(numElements, 0);
    blockSizes[0] = blockContents.size();
}

std::size_t Partition::getNumberOfBlocks() const {
    STORM_LOG_ASSERT(
        [this]() {
            uint64_t blocksCount = 0;
            forEachBlock([&blocksCount](const Block& block) { ++blocksCount; });
            return blocksCount;
        }() == numBlocks,
        "The cached number of blocks is inconsistent with the actual number of blocks");
    return numBlocks;
}

std::size_t Partition::getNumberOfElements() const {
    return elementToBlockIndex.size();
}

Partition::Block Partition::getUniversalBlock() const {
    return Block(blockContents.begin(), blockContents.end());
}

Partition::Block Partition::getBlockOfElement(ElementIndex element) const {
    STORM_LOG_ASSERT(element < elementToBlockIndex.size(), "Element index out of bounds");
    return getBlockFromIndex(elementToBlockIndex[element]);
}

bool Partition::contains(ElementIndex element, Block const& block) const {
    // We use the fact that each pair of blocks is either disjoint or one is a subset of the other.
    return isSubBlockOf(getBlockOfElement(element), block);
}

bool Partition::isSubBlockOf(Block const& subblock, Block const& superblock) const {
    STORM_LOG_ASSERT(checkBlockValidity(subblock), "Invalid subblock.");
    STORM_LOG_ASSERT(checkBlockValidity(superblock), "Invalid superblock.");
    // Check if subblock is a subspan of superblock
    return subblock.data() >= superblock.data() && (subblock.data() + subblock.size()) <= (superblock.data() + superblock.size());
}

bool Partition::isBlockOfElement(Block const& block, ElementIndex const& element) const {
    auto const& eBlock = getBlockOfElement(element);
    return eBlock.data() == block.data() && eBlock.size() == block.size();
}

Partition::BlockIndex Partition::getBlockIndex(Block const& block) const {
    STORM_LOG_ASSERT(checkBlockValidity(block), "Tried to get the index of an invalid block.");
    return std::distance(blockContents.data(), block.data());
}

Partition::Block Partition::getBlockFromIndex(BlockIndex blockIndex) const {
    STORM_LOG_ASSERT(blockIndex < blockContents.size(), "Block index out of bounds");
    STORM_LOG_ASSERT(blockSizes[blockIndex] > 0, "Block index points to a non-existent block");
    return Block(blockContents.data() + blockIndex, blockSizes[blockIndex]);
}

Partition::Block Partition::registerNewBlock(BlockIndex const start, BlockIndex const end) {
    STORM_LOG_ASSERT(start < end && end <= blockContents.size(), "Invalid block range");
    Block newBlock(blockContents.data() + start, end - start);
    // Shrinking an existing block doesn't require iterating over all block contents
    if (blockSizes[start] == 0) {
        for (ElementIndex const e : newBlock) {
            elementToBlockIndex[e] = start;
        }
    } else {
        STORM_LOG_ASSERT(isSubBlockOf(newBlock, getBlockFromIndex(start)), "New block is not a sub-block of the existing block");
    }
    blockSizes[start] = newBlock.size();
    return newBlock;
}

bool Partition::checkBlockValidity(Block const& block) const {
    // block must not be empty
    if (block.empty())
        return false;
    // the block refers to this partition
    if (block.data() < blockContents.data() || block.data() + block.size() > blockContents.data() + blockContents.size())
        return false;

    auto const blockIndex = std::distance(blockContents.data(), block.data());
    auto const blockEndIndex = blockIndex + block.size();

    // the block has the expected size and is either the last block or ends where another block starts
    if (!(blockSizes[blockIndex] == block.size() && (blockEndIndex == blockContents.size() || blockSizes[blockEndIndex] != 0)))
        return false;

    // The index of the first element must always match this block index (even if the element belongs to a subblock, it must have the same index)
    if (elementToBlockIndex[block[0]] != blockIndex)
        return false;

    // Check if the inverse mapping is consistent for all elements
    if (std::any_of(block.begin(), block.end(), [this](ElementIndex const& e) { return &blockContents[blockContentsInverse[e]] != &e; }))
        return false;

    if (isBlockOfElement(block, block[0])) {
        // this block has no proper sub-block so all elements should point to this blockIndex
        return std::all_of(block.begin(), block.end(),
                           [this, blockIndex, blockEndIndex](ElementIndex const& e) { return elementToBlockIndex[e] == blockIndex; });
    } else {
        // this block has proper sub-blocks. The block index of all elements should be within the range of this block.
        if (std::any_of(block.begin(), block.end(), [this, blockIndex, blockEndIndex](ElementIndex const& e) {
                // e's block index should be in the range of this block
                return !(elementToBlockIndex[e] >= blockIndex && elementToBlockIndex[e] < blockEndIndex);
            })) {
            return false;
        }
        // validate all sub-blocks
        bool allValid = true;
        forEachSubBlock(block, [this, &allValid](Block const& subBlock) { allValid = allValid && checkBlockValidity(subBlock); });
        return allValid;
    }
}

bool Partition::isProperSuperBlock(Block const& block) const {
    STORM_LOG_ASSERT(checkBlockValidity(block), "Tried to check if an invalid block is a proper superblock.");
    // block is a proper super-block iff it has a proper sub-block iff the block of some element is not equal to block
    return !isBlockOfElement(block, block[0]);
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

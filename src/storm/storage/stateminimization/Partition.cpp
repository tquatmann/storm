#include "storm/storage/stateminimization/Partition.h"

namespace storm::storage::stateminimization {

Partition::Partition(ElementIndex numElements) : numBlocks(1) {
    auto indexRange = std::ranges::iota_view<ElementIndex, ElementIndex>(0, numElements);
    blockContents.assign(indexRange.begin(), indexRange.end());
    blockContentsInverse = blockContents;
    elementToBlockIndex.assign(numElements, 0);
    blocks.resize(numElements);
    blocks[0] = Block(blockContents.data(), blockContents.size());
}

std::size_t Partition::getNumberOfBlocks() const {
    STORM_LOG_ASSERT([this]() {
        uint64_t blocksCount = 0;
        forEachBlock([&blocksCount](const Block& block) {
           ++blocksCount;
        });
        return blocksCount;
    }() == numBlocks, "The cached number of blocks is inconsistent with the actual number of blocks");
    return numBlocks;
}

std::size_t Partition::getNumberOfElements() const {
    return elementToBlockIndex.size();
}

Partition::Block const& Partition::getBlockOfElement(ElementIndex element) const {
    STORM_LOG_ASSERT(element < elementToBlockIndex.size(), "Element index out of bounds");
    return blocks[elementToBlockIndex[element]];
}


bool Partition::contains(ElementIndex element, Block const& block) const {
    // We use the fact that each pair of blocks is either disjoint or one is a subset of the other.
    return isSubBlockOf(getBlockOfElement(element), block);
}

bool Partition::isSubBlockOf(Block const& subblock, Block const& superblock) const {
    // Check if subblock is a subspan of superblock
    return subblock.data() >= superblock.data() &&
           (subblock.data() + subblock.size()) <= (superblock.data() + superblock.size());
}

bool Partition::isBlockOfElement(Block const& block, ElementIndex const& element) const {
    auto const& eBlock = getBlockOfElement(element);
    return eBlock.data() == block.data() && eBlock.size() == block.size();
}

Partition::BlockIndex Partition::getBlockIndex(Block const& block) const {
    STORM_LOG_ASSERT(checkBlockValidity(block), "Tried to get the index of an invalid block.");
    return std::distance(blockContents.data(), block.data());
}

// void Partition::invalidateCache(Partition::BlockIndex index) {
//     if (index < blockEndCache.size()) {
//         blockEndCache[index] = std::numeric_limits<BlockIndex>::max();
//     }
// }

bool Partition::checkBlockValidity(Block const& block) const {
    if (block.empty())
        return false;
    if (block.data() < blockContents.data() || block.data() + block.size() > blockContents.data() + blockContents.size())
        return false;

    // auto const blockIndex = std::distance(blockContents.data(), block.data());
    // auto const blockEndIndex = blockIndex + block.size();
    //
    // if (!blockIndices.get(blockIndex) || (blockEndIndex < blockIndices.size() && !blockIndices.get(blockEndIndex)))
    //     return false;
    // if (elementToBlockIndex[block[0]] != blockIndex)
    //     return false;
    //
    // return std::all_of(block.begin(), block.end(), [this, blockIndex, blockEndIndex](ElementIndex const& e) {
    //     // the given block might be a superblock, so e's block index does not need to coincide with blockIndex
    //     bool const isValidBlockIndex = elementToBlockIndex[e] >= blockIndex && elementToBlockIndex[e] <= blockEndIndex;
    //     return isValidBlockIndex  && blockContents[blockContentsInverse[e]] == e;
    // });
    return true;
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

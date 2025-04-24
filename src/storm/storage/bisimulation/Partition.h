#pragma once

#include <concepts>
#include <cstdint>
#include <iostream>
#include <set>
#include <span>
#include <type_traits>
#include <vector>

#include "storm/storage/BitVector.h"
#include "storm/utility/macros.h"

namespace storm::storage::bisimulation {

/*!
 * Represents a partition of a set of consecutive indices.
 * Allows efficient access to the blocks of the partition and the elements in each block.
 * Memory is cache-friendly and linear in the number of elements (independent of block count).
 * The downside is that merging blocks is not supported.
 */
class Partition {
   public:
    // Typedefs for readability
    using ElementIndex = uint64_t;
    using Block = std::span<ElementIndex const>;
    static constexpr auto BlockCompare = [](Block const& lhs, Block const& rhs) {
        if (lhs.data() < rhs.data()) {
            return true;
        }
        if (lhs.data() > rhs.data()) {
            return false;
        }
        return lhs.size() < rhs.size();
    };
    using BlockSet = std::set<Block, decltype(BlockCompare)>;

    Partition(Partition const& other) = default;
    Partition& operator=(Partition const& other) = default;
    Partition(Partition&& other) noexcept = default;
    Partition& operator=(Partition&& other) = default;

    /*!
     * Creates the partition { {0,1,...,numElements-1} } with a single block that contains all elements
     *
     * @param numElements the number of elements.
     */
    explicit Partition(ElementIndex numElements);

    /*!
     * Retrieves the number of blocks. Only counts those blocks that do not contain other (smaller) blocks
     */
    std::size_t getNumberOfBlocks() const;

    /*!
     * Retrieves the number of elements.
     */
    std::size_t getNumberOfElements() const;

    /*!
     * @return the block that contains the given element.
     * @note the lifetime of the returned object is limited to the lifetime of this partition.
     * @note subsequent operations on this partition (e.g. split) may change the order of the elements in the block, but do not affect the contents.
     * This means that changing the partition while iterating over block contents, i.e., `for (auto e : block) { partition.split(...) ... }` is not safe.
     */
    Block getBlockOfElement(ElementIndex element) const;

    /*!
     * @return true iff the given block contains the given element.
     * @note this has constant runtime, i.e., is preferrable over a linear search over the block
     */
    bool contains(Block const& block, ElementIndex element) const;

    /*!
     * Checks if the given block is a valid block in the partition. Usefull for sanity checks.
     * Will print an error if any of the checks fail, but does not assert/throw
     */
    bool checkBlockValidity(Block const& block) const;

    /*!
     * @return true iff the given block contains multiple sub-block. For example, {1,2,3} is a (proper) super block of {1,2} and {3}.
     */
    bool isProperSuperBlock(Block const& block) const;

    /*!
     * Applies the given function on each block in the partition.
     * @note: Splitting the currently processed block in the function *is* valid. However, the resulting sub-blocks are not processed recursively.
     */
    template<typename Func>
        requires std::invocable<Func, Block>
    void forEachBlock(Func const& f) const {
        forEachBlockInIndexRange(0u, blockIndices.size(), f);
    }

    /*!
     * Applies the given function on each block in the given super block.
     * For example, if superBlock = {0,1,2,3,4} and the partition currently has blocks {{0}, {1,2,3}, {4}, ...},
     * the function will be called with the blocks {0}, {1,2,3}, and {4}.
     * @note: Splitting the currently processed block in the function *is* valid. However, the resulting sub-blocks are not processed recursively.
     */
    template<typename Func>
        requires std::invocable<Func, Block>
    void forEachSubBlock(Block const& superBlock, Func const& f) const {
        STORM_LOG_ASSERT(!superBlock.empty(), "Superblock is empty.");
        auto const blockStart = getBlockIndex(superBlock);
        auto const blockEnd = blockStart + superBlock.size();
        forEachBlockInIndexRange(blockStart, blockEnd, f);
    }

    /*!
     * Splits the given block according to the given order.
     * Specifically, the elements in the block are sorted according to the given order and then the block is split into
     * multiple blocks, divided at every position where the order changes.
     * @return true iff the block was split, i.e. if the input block is now a proper super block.
     */
    template<typename SplittingOrder>
        requires std::invocable<SplittingOrder, ElementIndex, ElementIndex>
    bool splitBlockByOrder(Block const& block, SplittingOrder const& less) {
        if (block.size() <= 1) {
            return false;  // nothing to do
        }

        // Sort the contents of the block first
        STORM_LOG_ASSERT(!isProperSuperBlock(block), "Tried to split a block that consists of multiple sub-blocks.");
        auto const blockStart = getBlockIndex(block);
        auto const blockEnd = blockStart + block.size();
        std::sort(blockContents.begin() + blockStart, blockContents.begin() + blockEnd, less);

        // Catch the special case where there is no split
        if (!less(blockContents[blockStart], blockContents[blockEnd - 1])) {
            return false;  // nothing to do
        }

        // helper function to find the end index of a current block
        auto getEndOfBlock = [this, less, blockEnd](BlockIndex const currIndex) {
            for (auto i = currIndex + 1; i < blockEnd; ++i) {
                if (less(blockContents[currIndex], blockContents[i])) {
                    return i;
                }
            }
            return blockEnd;
        };

        // Now create new blocks whenever the order changes
        for (auto newBlockIndex = getEndOfBlock(blockStart); newBlockIndex < blockEnd;) {
            STORM_LOG_ASSERT(!blockIndices.get(newBlockIndex),
                             "Partition in inconsistent state: Block index " << newBlockIndex << " already set as start index.");
            blockIndices.set(newBlockIndex);
            auto const newBlockEnd = getEndOfBlock(newBlockIndex);
            std::for_each(blockContents.begin() + newBlockIndex, blockContents.begin() + newBlockEnd,
                          [this, &newBlockIndex](ElementIndex const& e) { elementToBlockIndex[e] = newBlockIndex; });
            newBlockIndex = newBlockEnd;
        }
        return true;  // there must have been a split because the case without a split is already  catched above
    }

    /*!
     * Splits the given block according to the given predicate.
     * Specifically, the elements in the block are swapped around, such that all elements where the predicate f evaluates to false come first.
     * Then, two sub-blocks are created accordingly. If all elements are true/false, no splitting is performed.
     * @return a pair containing first the 'false' sub-block and then the 'true' sub-block. One of them can be empty.
     */
    template<typename SplittingPredicate>
        requires std::predicate<SplittingPredicate, ElementIndex>
    std::pair<Block, Block> splitBlockByPredicate(Block const& block, SplittingPredicate const& f) {
        STORM_LOG_ASSERT(!isProperSuperBlock(block), "Tried to split a block that consists of multiple sub-blocks.");
        STORM_LOG_ASSERT(!block.empty(), "Tried to split an empty block");

        // swap the block contents so that all elements e with f(e) == false come first
        auto const blockStart = getBlockIndex(block);
        auto const blockEnd = blockStart + block.size();
        auto l = blockStart;
        auto r = blockEnd - 1;
        // Loop invariant: all elements at position < l are false, all elements at position > r are true
        while (l <= r) {
            while (l <= r && !f(std::as_const(blockContents[l]))) {
                ++l;
            }
            if (l > r) {
                break;  // no more elements to swap
            }
            // At this point we know that f(blockContents[l]) == true
            while (l < r && f(std::as_const(blockContents[r]))) {
                --r;
            }
            if (l == r) {
                // We have f(blockContents[r]) == f(blockContents[l]) == true
                --r;
                break;  // l > r holds now
            } else if (l < r) {
                std::swap(blockContents[l], blockContents[r]);
                ++l;
                --r;
            }
        }
        STORM_LOG_ASSERT(l == r + 1, "Unexpected indices");

        // Handle cases where there is no split
        if (l == blockStart) {
            return {Block{}, block};  // all elements are true
        } else if (l == blockEnd) {
            return {block, Block{}};  // all elements are false
        }

        // Perform a split
        auto const newBlockIndex = l;
        blockIndices.set(newBlockIndex);
        std::for_each(blockContents.begin() + newBlockIndex, blockContents.begin() + blockEnd,
                      [this, &newBlockIndex](ElementIndex const& e) { elementToBlockIndex[e] = newBlockIndex; });
        return {getBlockFromIndexRange(blockStart, newBlockIndex), getBlockFromIndexRange(newBlockIndex, blockEnd)};
    }

   private:
    /*!
     * The index of a block in the partition.
     * @note we make this private since we do not want to expose the internal representation of the partition.
     */
    using BlockIndex = uint64_t;

    /*!
     * @return the block indices of the given block.
     * @note block indices are invalidated when the partition is modified (e.g. after calling split).
     * We therefore do not want to expose them.
     */
    BlockIndex getBlockIndex(Block const& block) const;

    /*!
     * Creates a block from the given range of indices.
     * @return a (super?)block that contains all states whose block index is in the range [start, end).
     */
    Block getBlockFromIndexRange(BlockIndex const start, BlockIndex const end) const;

    /*!
     * Creates a block from the given range of indices.
     * @return a (super?)block that contains all states whose block index is in the range [start, end).
     */
    Block getBlockFromIndex(BlockIndex const index) const;

    /*!
     * Iterates over all blocks in the given range of block indices [start, end) and applies the given function to each block.
     * @note we make this private since we do not want to expose the internal representation of the partition.
     * @note: Splitting the currently processed block in the function *is* valid. However, the resulting sub-blocks are not processed recursively.
     */
    template<typename Func>
        requires std::invocable<Func, Block>
    void forEachBlockInIndexRange(BlockIndex const start, BlockIndex const end, Func const& f) const {
        for (auto blockIndex = start; blockIndex < end;) {
            auto block = getBlockFromIndex(blockIndex);
            blockIndex += block.size();
            f(block);
        }
    }

    /// Stores for each block the elements in that block (cf. blockIndices)
    /// Stores where a new block begins in the blockContents vector.
    /// The number of set bits equals the number of blocks. The first bit is always set.
    /// The k'th block starts at the k'th set bit. The BlockIndex of the k'th block is the position of that bit.
    /// If bit i is set, the corresponding block is given by { blockContents[j] | i ≤ j < blockIndices.getNextSetIndex(i+1) }
    std::vector<ElementIndex> blockContents;
    storm::storage::BitVector blockIndices;

    /// Maps each element to the start index of its block.
    /// for all elements s, blockIndices.get(elementToBlockIndex[s]) is true and s is in { blockContents[j] | elementToBlockIndex[s] ≤ j <
    /// blockIndices.getNextSetIndex(elementToBlockIndex[s]+1) }
    std::vector<BlockIndex> elementToBlockIndex;
};

std::ostream& operator<<(std::ostream& os, const Partition& partition);

}  // namespace storm::storage::bisimulation

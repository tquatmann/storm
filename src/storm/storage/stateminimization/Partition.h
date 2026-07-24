#pragma once

#include <boost/sort/block_indirect_sort/block_indirect_sort.hpp>
#include <cstddef>
#include <memory>
#include <ranges>

#include "storm/storage/BitVector.h"
#include "storm/utility/macros.h"

namespace storm::storage::stateminimization {

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

    struct BlockCompare {
        bool operator()(Block const& lhs, Block const& rhs) const {
            if (lhs.size() < rhs.size()) {
                return true;
            }
            if (lhs.size() > rhs.size()) {
                return false;
            }
            return lhs.data() < rhs.data();
        }
    };
    // Arbitrary set of blocks, ordered by their size (smallest first)
    using OrderedBlockSet = std::set<Block, BlockCompare>;

    /*!
     * Set of blocks that are no proper superblock.
     * It is illegal to split a block that is currently in this set.
     * Insertion and popping an arbitrary element has mostly constant time.
     * However, the set requires Theta(numElements) space.
     */
    class NonSuperBlockSet {
        public:
        NonSuperBlockSet(Partition const& partition) : partition(partition), blockIndices(partition.getNumberOfElements(), false) {}

        std::size_t size() const {
            return containedBlockIndices.size();
        }
        bool empty() const {
            return containedBlockIndices.empty();
        }
        bool contains(Block const& block) const {
            return blockIndices.get(partition.getBlockIndex(block));
        }

        void insert(Block const& block) {
            STORM_LOG_ASSERT(!block.empty(), "Cannot insert empty block.");
            STORM_LOG_ASSERT(!partition.isProperSuperBlock(block), "Cannot insert block that is a proper superblock.");
            if (BlockIndex i = partition.getBlockIndex(block); !blockIndices.get(i)) {
                blockIndices.set(i, true);
                containedBlockIndices.push_back(i);
            }
        }

        Block pop() {
            STORM_LOG_ASSERT(!empty(), "Cannot pop from empty blockset.");
            auto const i = containedBlockIndices.back();
            blockIndices.set(i, false);
            containedBlockIndices.pop_back();
            return partition.getBlockFromIndex(i);
        }



    private:

        Partition const& partition;
        std::vector<uint64_t> containedBlockIndices;
        storm::storage::BitVector blockIndices;
    };

    Partition() = default;
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
     * @return the block that contains all elements, i.e., the universal block.
     */
Block getUniversalBlock() const {
        return Block(blockContents.begin(), blockContents.end());
    }

    /*!
     * @return the block that contains the given element.
     * @note the lifetime of the returned object is limited to the lifetime of this partition.
     * @note subsequent operations on this partition (e.g. split) may change the order of the elements in the block, but do not affect the contents.
     * This means that changing the partition while iterating over block contents, i.e., `for (auto e : block) { partition.split(...) ... }` is not safe.
     */
    Block getBlockOfElement(ElementIndex element) const;

    /*!
     * @return true iff the given element is contained in the given block
     * @note this has constant runtime, i.e., is preferrable over a linear search over the block
     */
    bool contains(ElementIndex element, Block const& block) const;

    /*!
     * @return true iff all elements of the given subblock are contained in the given superblock
     * @note this has constant runtime, i.e., is preferrable over a linear search over the block
     */
    bool isSubBlockOf(Block const& subblock, Block const& superblock) const;


    /*!
 * @return true iff the (smallest known) block of the given element coincides with the given block
 */
    bool isBlockOfElement(Block const& block, ElementIndex const& element) const;


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
        forEachSubBlock(getUniversalBlock(), f);
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
        auto superBlockEnd = getBlockIndex(superBlock) + superBlock.size();
        for (BlockIndex blockStart = getBlockIndex(superBlock); blockStart < superBlockEnd;) {
            auto const& subBlock = getBlockOfElement(blockContents[blockStart]);
            blockStart += subBlock.size();
            f(subBlock);
        };
    }

    template<typename Iterator, typename Compare>
    void small_sort(Iterator begin, Iterator end, Compare comp) {
        // STORM_LOG_ASSERT(false, "try other sorting"); TODO
        if (std::distance(begin, end) <= 32) {
            for (auto i = begin + 1; i != end; ++i) {
                auto key = *i;
                auto j = i;
                while (j > begin && comp(key, *(j - 1))) {
                    *j = *(j - 1);
                    --j;
                }
                *j = key;
            }
        } else {
            boost::sort::block_indirect_sort(begin, end, comp);
        }
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
        // small_sort(blockContents.begin() + blockStart, blockContents.begin() + blockEnd, less);
        std::sort(blockContents.begin() + blockStart, blockContents.begin() + blockEnd, less);

        // update the inverse after sorting
        for (ElementIndex i = blockStart; i < blockEnd; ++i) {
            blockContentsInverse[blockContents[i]] = i;
        }

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
        for (auto newBlockIndex = blockStart; newBlockIndex < blockEnd;) {
            auto const newBlockEnd = getEndOfBlock(newBlockIndex);
            registerNewBlock(newBlockIndex, newBlockEnd);
            newBlockIndex = newBlockEnd;
            ++numBlocks;
        }
        --numBlocks; // the old superblock is no longer a block as counted by this partition
        STORM_LOG_ASSERT(isProperSuperBlock(block), "Partition in inconsistent state: Block was not split into multiple sub-blocks.");
        return true;  // there must have been a split because the case without a split is already caught above
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
                blockContentsInverse[blockContents[l]] = l;
                blockContentsInverse[blockContents[r]] = r;
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
        Block noBlock = registerNewBlock(blockStart, l);
        Block yesBlock = registerNewBlock(l, blockEnd);
        ++numBlocks; // We have split a single block into two blocks, so the total number increments by 1
        return {noBlock, yesBlock};
    }

    /*!
     * Splits the given block according to the given range of elements.
     * Specifically, the elements in the block are swapped around, such that all elements that are not in the given range come first.
     * Then, two sub-blocks are created accordingly. If the range is empty or contains all elements of the block, no splitting is performed.
     * @param r the input range.
     * @return a pair containing first the sub-block that is not in the range and then the sub-block in the range. One of them can be empty.
     */
    template<typename SplitRange>
        requires std::ranges::input_range<SplitRange> && std::same_as<std::ranges::range_value_t<SplitRange>, ElementIndex>
    std::pair<Block, Block> splitBlockByRange(Block const& block, SplitRange const& r) {
        STORM_LOG_ASSERT(!isProperSuperBlock(block), "Tried to split a block that consists of multiple sub-blocks.");
        STORM_LOG_ASSERT(!block.empty(), "Tried to split an empty block");

        auto const blockStart = getBlockIndex(block);
        auto const blockEnd = blockStart + block.size();

        // We split the block into a "no"-part and a "yes"-part

        auto yesBlockStart = blockEnd;
        for (auto const& element : r) {
            auto const src = blockContentsInverse[element];
            if (src < blockStart || src >= yesBlockStart) {
                continue; // element is either not in the block or already occurred before (duplicate in r)
            }
            // swap the element to the end of the block, into the 'yes' part
            --yesBlockStart;
            std::swap(blockContents[src], blockContents[yesBlockStart]);
            blockContentsInverse[blockContents[src]] = src;
            blockContentsInverse[blockContents[yesBlockStart]] = yesBlockStart;
        }

        // Handle cases where there is no split
        if (yesBlockStart == blockStart) {
            return {Block{}, block};  // all elements are in the range
        } else if (yesBlockStart == blockEnd) {
            return {block, Block{}};  // no elements are in the range
        }

        // Perform a split
        Block noBlock = registerNewBlock(blockStart, yesBlockStart);
        Block yesBlock = registerNewBlock(yesBlockStart, blockEnd);
        ++numBlocks; // We have split a single block into two blocks, so the total number increments by 1
        return {noBlock, yesBlock};
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

    Block getBlockFromIndex(BlockIndex blockIndex) const {
        STORM_LOG_ASSERT(blockIndex < blockContents.size(), "Block index out of bounds");
        STORM_LOG_ASSERT(blockSizes[blockIndex] > 0, "Block index points to a non-existend block");
        return Block(blockContents.data() + blockIndex, blockSizes[blockIndex]);
    }

    /*!
     * Creates a block from the given range of indices.
     * @return a (super?)block that contains all states whose block index is in the range [start, end).
     */

    Block registerNewBlock(BlockIndex const start, BlockIndex const end) {
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


    /// Stores for each block the elements in that block (cf. blockIndices)
    /// Stores where a new block begins in the blockContents vector.
    /// The number of set bits equals the number of blocks. The first bit is always set.
    /// The k'th block starts at the k'th set bit. The BlockIndex of the k'th block is the position of that bit.
    /// If bit i is set, the corresponding block is given by { blockContents[j] | i ≤ j < blockIndices.getNextSetIndex(i+1) }
    std::vector<ElementIndex> blockContents;
    std::vector<BlockIndex> blockContentsInverse;

    /// Maps each element to the start index of its block.
    /// for all elements s, blockIndices.get(elementToBlockIndex[s]) is true and s is in { blockContents[j] | elementToBlockIndex[s] ≤ j <
    /// blockIndices.getNextSetIndex(elementToBlockIndex[s]+1) }
    std::vector<BlockIndex> elementToBlockIndex;
    std::vector<std::size_t> blockSizes;

    uint64_t numBlocks;

};

std::ostream& operator<<(std::ostream& os, const Partition& partition);

}  // namespace storm::storage::stateminimization

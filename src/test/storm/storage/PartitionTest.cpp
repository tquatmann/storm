#include <sstream>
#include "storm/storage/bisimulation/Partition.h"

#include "test/storm_gtest.h"

namespace {

std::string blockToString(auto const& block) {
    std::ostringstream oss;
    oss << "{";
    for (bool first{true}; auto const e : block) {
        if (first) {
            first = false;
        } else {
            oss << ", ";
        }
        oss << e;
    }
    oss << "}";
    return oss.str();
}

bool equalBlocks(storm::storage::bisimulation::Partition const& partition, std::initializer_list<uint64_t> const& expected,
                 storm::storage::bisimulation::Partition::Block const& actual) {
    bool const expectedContainsActual =
        std::all_of(actual.begin(), actual.end(), [&expected](auto const& e) { return std::find(expected.begin(), expected.end(), e) != expected.end(); });
    EXPECT_TRUE(expectedContainsActual) << "Content missmatch: " << blockToString(expected) << " vs. " << blockToString(actual);
    bool const actualContainsExpected =
        std::all_of(expected.begin(), expected.end(), [&actual](auto const& e) { return std::find(actual.begin(), actual.end(), e) != actual.end(); });
    EXPECT_TRUE(actualContainsExpected) << "Content missmatch: " << blockToString(expected) << " vs. " << blockToString(actual);

    // Also test the contains method
    for (auto const& e : expected) {
        EXPECT_TRUE(partition.contains(actual, e)) << "Element " << e << " not in block " << blockToString(actual);
    }

    return expectedContainsActual && actualContainsExpected;
}

TEST(PartitionTest, Basic) {
    storm::storage::bisimulation::Partition partition(10);

    // Check initial partition
    EXPECT_EQ(10ul, partition.getNumberOfElements());
    EXPECT_EQ(1ul, partition.getNumberOfBlocks());
    for (auto i = 0ul; i < 10; ++i) {
        auto block = partition.getBlockOfElement(i);
        ASSERT_TRUE(partition.checkBlockValidity(block)) << "Block " << blockToString(block) << " is not valid.";
        EXPECT_TRUE(equalBlocks(partition, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, block));
    }
    auto initBlock = partition.getBlockOfElement(0);
    ASSERT_FALSE(partition.isProperSuperBlock(initBlock));

    // Split the partition into two blocks (odd and even)
    partition.forEachBlock([&partition](auto const& block) {
        auto [odd, even] = partition.splitBlockByPredicate(block, [](auto const& e) { return e % 2 == 0; });
        ASSERT_TRUE(partition.checkBlockValidity(odd)) << "Block " << blockToString(odd) << " is not valid.";
        ASSERT_TRUE(partition.checkBlockValidity(even)) << "Block " << blockToString(even) << " is not valid.";
        EXPECT_TRUE(equalBlocks(partition, {1, 3, 5, 7, 9}, odd));
        EXPECT_TRUE(equalBlocks(partition, {0, 2, 4, 6, 8}, even));
    });
    ASSERT_TRUE(partition.checkBlockValidity(initBlock)) << "Block " << blockToString(initBlock) << " is not valid.";
    EXPECT_EQ(2ul, partition.getNumberOfBlocks());
    ASSERT_TRUE(partition.isProperSuperBlock(initBlock));
    auto oddBlock = partition.getBlockOfElement(1);
    auto evenBlock = partition.getBlockOfElement(0);

    // Split the odd block w.r.t. their  remainder when divided by 3
    EXPECT_TRUE(partition.splitBlockByOrder(oddBlock, [](auto const& a, auto const& b) { return a % 3 < b % 3; }));
    uint64_t seenSubBlocks = 0;  // very "smart" (aka lazy) encoding of the seen subsets
    partition.forEachSubBlock(oddBlock, [&partition, &seenSubBlocks](auto const& block) {
        ASSERT_TRUE(partition.checkBlockValidity(block)) << "Block " << blockToString(block) << " is not valid.";
        if (partition.contains(block, 1)) {
            EXPECT_TRUE(equalBlocks(partition, {1, 7}, block));
            seenSubBlocks += 1;
        } else if (partition.contains(block, 3)) {
            EXPECT_TRUE(equalBlocks(partition, {3, 9}, block));
            seenSubBlocks += 10;
        } else if (partition.contains(block, 5)) {
            EXPECT_TRUE(equalBlocks(partition, {5}, block));
            seenSubBlocks += 100;
        } else {
            FAIL() << "Block " << blockToString(block) << " does not contain any expected elements.";
        }
    });
    EXPECT_EQ(111ul, seenSubBlocks) << "Not all sub blocks have been seen.";
    EXPECT_EQ(4ul, partition.getNumberOfBlocks());

    // Check if the contains method works as expected
    for (auto i = 0ul; i < 10; ++i) {
        EXPECT_EQ(partition.contains(evenBlock, i), i % 2 == 0) << "Element " << i << " not in block " << blockToString(i % 2 == 0 ? evenBlock : oddBlock);
        EXPECT_EQ(partition.contains(oddBlock, i), i % 2 != 0) << "Element " << i << " in block " << blockToString(i % 2 != 0 ? oddBlock : evenBlock);
    }

    EXPECT_TRUE(partition.isProperSuperBlock(oddBlock));
    EXPECT_TRUE(equalBlocks(partition, {1, 3, 5, 7, 9}, oddBlock));

    // Perform splits that have no effect, i.e., do not change the number of blocks
    auto [evenAndOddBlock, evenAndEvenBlock] = partition.splitBlockByPredicate(evenBlock, [](auto const& e) { return e % 2 == 0; });
    EXPECT_TRUE(evenAndOddBlock.empty());
    EXPECT_TRUE(equalBlocks(partition, {0, 2, 4, 6, 8}, evenAndEvenBlock));
    EXPECT_FALSE(partition.isProperSuperBlock(evenBlock));

    auto [less10, geq10] = partition.splitBlockByPredicate(evenBlock, [](auto const& e) { return e >= 10; });
    EXPECT_TRUE(geq10.empty());
    EXPECT_TRUE(equalBlocks(partition, {0, 2, 4, 6, 8}, less10));
    EXPECT_FALSE(partition.isProperSuperBlock(evenBlock));

    EXPECT_FALSE(partition.splitBlockByOrder(evenBlock, [](auto const& a, auto const& b) { return a % 2 < b % 2; }));
    EXPECT_FALSE(partition.isProperSuperBlock(evenBlock));

    typename storm::storage::bisimulation::Partition::BlockSet blockSet;
    blockSet.insert(evenBlock);
    blockSet.insert(oddBlock);
    EXPECT_EQ(blockSet.size(), 2);
    blockSet.insert(evenAndEvenBlock);  // (re-)inserts the even block
    EXPECT_EQ(blockSet.size(), 2);
    blockSet.insert(partition.getBlockOfElement(0));  // (re-)inserts the even block
    EXPECT_EQ(blockSet.size(), 2);
    EXPECT_EQ(blockSet.size(), 2);
    blockSet.insert(partition.getBlockOfElement(1));  // inserts a sub-block of oddBlock
    EXPECT_EQ(blockSet.size(), 3);
};
}  // namespace

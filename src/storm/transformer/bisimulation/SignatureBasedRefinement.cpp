#include "storm/transformer/bisimulation/SignatureBasedRefinement.h"

#include <algorithm>
#include <cstdint>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Signatures.h"
#include "storm/transformer/bisimulation/SparseAccumulator.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

namespace detail {

template<typename ValueType, SignatureMode SignatureMode, typename QuotientValueType>
struct SignatureRefinementContext {
    SignatureRefinementContext(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                               Signatures<ValueType, SignatureMode, QuotientValueType>& signatures)
        : model(model), partition(partition), signatures(signatures), backwardTransitions(model.getBackwardTransitions()), cache(partition) {}

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition& partition;
    storm::bisimulation::Signatures<ValueType, SignatureMode, QuotientValueType>& signatures;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    storm::bisimulation::Partition::OrderedBlockMap<bool>
        queue;  // stores an extra flag for each element in the queue. The flag indicates whether we enforce exploring the predecessors of the block

    struct Cache {
        Cache(Partition const& partition) : predecessorToPivotBlocks(partition.getNumberOfElements()), predecessorBlocks(partition) {}
        WrappingSetAccumulator predecessorToPivotBlocks;
        Partition::NonSuperBlockSet predecessorBlocks;
    } cache;
};

template<typename ValueType, SignatureMode SignatureMode, typename QuotientValueType>
void refinePartitionBasedOnSignature(SignatureRefinementContext<ValueType, SignatureMode, QuotientValueType>& context,
                                     storm::bisimulation::Partition::Block const pivotBlock, bool const enforcePredecessorExploration) {
    // Split the pivot block B into B=B_1 cup B_2 cup ... cup B_n using signature refinement
    // First update the state signatures
    for (uint64_t const state : pivotBlock) {
        context.signatures.updateStateSignature(state);
    }
    // Then perform the signature-based split.
    bool pivotHasBeenSplit{false};
    if constexpr (SignatureMode == storm::bisimulation::SignatureMode::Exact) {
        pivotHasBeenSplit = context.partition.splitBlockByOrder(pivotBlock, context.signatures.getEquivalenceSplitOrder());
    } else {
        // In approximative mode or interval abstraction mode, there is no transitive order on signatures that captures "approximately equal signatures".
        // We therefore split in two steps: first, we split by a coarse order based on structural properties of the state signature. Then, we do a more
        // expensive split based on clustering on each sub-block.
        pivotHasBeenSplit = context.partition.splitBlockByOrder(pivotBlock, context.signatures.getStructuralSplitOrder());
        auto splitCondition = context.signatures.getClusteringSplitCondition();
        context.partition.forEachSubBlock(pivotBlock, [&context, &splitCondition, &pivotHasBeenSplit](auto const& subBlock) {
            pivotHasBeenSplit |= context.partition.splitBlockByClustering(subBlock, splitCondition);
        });
        // In interval-abstraction mode, the clustering-based split also triggers widening of the signature of the *first* state of each subBlock so that it contains all behaviour of the remaining state signatures in the subblock.
        // Once the partition is final, we will ensure that every state of a block gets that representative signature. See performSignatureBasedRefinement.
    }

    if (!pivotHasBeenSplit && !enforcePredecessorExploration) {
        // When the current pivot block is stable, there is no need to look into its predecessors. We can continue with the next pivot.
        return;
    }

    // When the pivot block is not stable, it means it has been split and the predecessor blocks have not been checked since that split.
    // Therefore, all the predecessor blocks need to checked again.
    // While we could just add all those predecessors to the queue, we instead try to split them first based on simple, graph-based criteria so that (expensive)
    // signature refinement is hopefully only applied to smaller blocks. Specifically, we split predecessor blocks based on which set of sub-blocks of the pivot
    // they can reach.

    // Gather predecessors of the pivot and their reachable sub-blocks
    auto& predecessorToPivotBlocks = context.cache.predecessorToPivotBlocks;
    auto& predecessorBlocks = context.cache.predecessorBlocks;
    uint64_t subBlockIndex = 0;
    context.partition.forEachSubBlock(pivotBlock, [&context, &predecessorToPivotBlocks, &predecessorBlocks, &subBlockIndex](auto const& subBlock) {
        for (uint64_t const state : subBlock) {
            for (auto const& predecessorEntry : context.backwardTransitions.getRow(state)) {
                auto const predecessorState = predecessorEntry.getColumn();
                auto const predecessorBlock = context.partition.getBlockOfElement(predecessorState);
                if (predecessorBlock.size() > 1) {
                    // No need to investigate singleton predecessor blocks as they cannot be split any further.
                    predecessorToPivotBlocks.addValue(predecessorState, subBlockIndex);
                    predecessorBlocks.insert(predecessorBlock);
                }
            }
        }
        ++subBlockIndex;
    });
    // Apply splitting of the pivot predecessors (similar to splitter-based refinement)
    while (!predecessorBlocks.empty()) {
        auto const predecessorBlock = predecessorBlocks.pop();
        // Split the predecessor block according to which sub-blocks of the pivot-block can be reached.
        auto const& toPivotBlocks = predecessorToPivotBlocks.getValues();

        // Split the block by whether a state is a predecessor of the pivot block or not
        // We do this by either iterating over the pivotPredecessorStates or the predecessorBlock, depending on what is shorter.
        auto [noPredecessors, predecessors] =
            predecessorToPivotBlocks.getNonDefaultStates().size() < predecessorBlock.size()
                ? context.partition.splitBlockByRange(predecessorBlock, predecessorToPivotBlocks.getNonDefaultStates())
                : context.partition.splitBlockByPredicate(predecessorBlock, [&toPivotBlocks](auto const& state) { return toPivotBlocks[state] != 0ull; });

        // At least one state should be a predecessor of the pivot block (otherwise we wouldn't have found that block above)
        STORM_LOG_ASSERT(!predecessors.empty(), "The predecessor block should contain at least one predecessor state.");

        // Now apply the splitting based on which sub-blocks of the pivot block can be reached.
        // If we did not actually split the pivot block, this operation would have no effect.
        if (pivotHasBeenSplit) {
            context.partition.splitBlockByOrder(
                predecessors, [&toPivotBlocks](uint64_t const state1, uint64_t const state2) { return toPivotBlocks[state1] < toPivotBlocks[state2]; });
        } else {
            STORM_LOG_ASSERT(
                std::all_of(predecessors.begin(), predecessors.end(), [&toPivotBlocks](uint64_t const& state) { return toPivotBlocks[state] == 1ull; }),
                "Expected all predecessor states to reach the pivot block.");
        }

        if (context.partition.isProperSuperBlock(predecessorBlock)) {
            // Erase the super block from the queue and the enforced unstable blocks (it might or might not be in there)
            context.queue.erase(predecessorBlock);
            // Add all sub-blocks to the queue. As we made a split, we must explore the predecessors of predecessorBlock.
            context.partition.forEachSubBlock(predecessorBlock, [&context](auto const& block) { context.queue[block] = true; });
        } else {
            // The simple, graph-based splitting was not effective. We must add the entire predecessorBlock to the queue. We do not have to
            // enforce that predecessors are explored.
            context.queue.try_emplace(predecessorBlock, false);
        }
    }
    // Reset the touched reachable-subblocks
    predecessorToPivotBlocks.clear();
}

}  // namespace detail

template<typename ValueType, SignatureMode SignatureMode, typename QuotientValueType>
void performSignatureBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                     Signatures<ValueType, SignatureMode, QuotientValueType>& signatures) {
    detail::SignatureRefinementContext<ValueType, SignatureMode, QuotientValueType> context(model, partition, signatures);
    // Initially, add all current blocks to the queue. No need to enforce exploring predecessors.
    partition.forEachBlock([&context](auto const& block) { context.queue.emplace(block, false); });

    while (!context.queue.empty()) {
        // take the smallest block from the queue
        auto const [pivotBlock, enforcePredecessorExploration] = *context.queue.begin();
        context.queue.erase(context.queue.begin());
        STORM_LOG_ASSERT(!partition.isProperSuperBlock(pivotBlock), "Broken invariant: the queue should not contain blocks that have been split.");
        // Split the pivotBlock based on its signature and split the predecessor blocks based on a simple, graph-based criterion
        detail::refinePartitionBasedOnSignature(context, pivotBlock, enforcePredecessorExploration);
    }

    // Singleton blocks are never re-examined once formed (because they cannot be split any further). Their cached signature can become stale if one of their
    // successor blocks is split afterwards. Refresh them here so that, as documented, all signatures are up to date once this function returns.
    partition.forEachBlock([&signatures](auto const& block) {
        if (block.size() == 1) {
            signatures.updateStateSignature(block.front());
        } else if constexpr (SignatureMode == bisimulation::SignatureMode::IntervalAbstraction) {
            // refinePartitionBasedOnSignature only ensures that the widened (abstracted) intervals are stored at the *first* state of each block. We copy it over to the other states, too.
            for (auto blockIt = block.begin() + 1; blockIt != block.end(); ++blockIt) {
                signatures.copyStructuralEquivalentStateSignature(block.front(), *blockIt);
            }
        }
    });
}

// Explicit instantiations for QuotientValueType == ValueType
template void performSignatureBasedRefinement<double, SignatureMode::Exact>(storm::models::sparse::Model<double> const& model,
                                                                            storm::bisimulation::Partition& partition,
                                                                            Signatures<double, SignatureMode::Exact>& signatures);
template void performSignatureBasedRefinement<double, SignatureMode::Approximative>(storm::models::sparse::Model<double> const& model,
                                                                                    storm::bisimulation::Partition& partition,
                                                                                    Signatures<double, SignatureMode::Approximative>& signatures);
template void performSignatureBasedRefinement<storm::RationalNumber, SignatureMode::Exact>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                                           storm::bisimulation::Partition& partition,
                                                                                           Signatures<storm::RationalNumber, SignatureMode::Exact>& signatures);
template void performSignatureBasedRefinement<storm::RationalNumber, SignatureMode::Approximative>(
    storm::models::sparse::Model<storm::RationalNumber> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::RationalNumber, SignatureMode::Approximative>& signatures);
template void performSignatureBasedRefinement<storm::RationalFunction, SignatureMode::Exact>(
    storm::models::sparse::Model<storm::RationalFunction> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::RationalFunction, SignatureMode::Exact>& signatures);
template void performSignatureBasedRefinement<storm::Interval, SignatureMode::Exact>(storm::models::sparse::Model<storm::Interval> const& model,
                                                                                     storm::bisimulation::Partition& partition,
                                                                                     Signatures<storm::Interval, SignatureMode::Exact>& signatures);
template void performSignatureBasedRefinement<storm::Interval, SignatureMode::IntervalAbstraction>(
    storm::models::sparse::Model<storm::Interval> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::Interval, SignatureMode::IntervalAbstraction>& signatures);
template void performSignatureBasedRefinement<storm::RationalInterval, SignatureMode::Exact>(
    storm::models::sparse::Model<storm::RationalInterval> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::RationalInterval, SignatureMode::Exact>& signatures);
template void performSignatureBasedRefinement<storm::RationalInterval, SignatureMode::IntervalAbstraction>(
    storm::models::sparse::Model<storm::RationalInterval> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::RationalInterval, SignatureMode::IntervalAbstraction>& signatures);

// Explicit instantiations for QuotientValueType == IntervalType<ValueType> (for interval abstraction)
// Exact mode is not meaningful in this case, as that would mean that we are never allowed to abstract values into intervals.
template void performSignatureBasedRefinement<double, SignatureMode::IntervalAbstraction, storm::Interval>(
    storm::models::sparse::Model<double> const& model, storm::bisimulation::Partition& partition,
    Signatures<double, SignatureMode::IntervalAbstraction, storm::Interval>& signatures);
template void performSignatureBasedRefinement<storm::RationalNumber, SignatureMode::IntervalAbstraction, storm::RationalInterval>(
    storm::models::sparse::Model<storm::RationalNumber> const& model, storm::bisimulation::Partition& partition,
    Signatures<storm::RationalNumber, SignatureMode::IntervalAbstraction, storm::RationalInterval>& signatures);

}  // namespace storm::bisimulation

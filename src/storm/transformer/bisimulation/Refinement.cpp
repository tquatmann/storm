#include "storm/transformer/bisimulation/Refinement.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

namespace detail {

template<typename ValueType>
struct RefinementContext {
    RefinementContext(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                      storm::IntervalBaseType<ValueType> tolerance)
        : model(model), partition(partition), backwardTransitions(model.getBackwardTransitions()), tolerance(tolerance) {}

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition& partition;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    ValueType const tolerance;
    storm::bisimulation::Partition::OrderedBlockSet queue;
};

template<typename ValueType>
struct SplitterRefinementCache {
    std::vector<ValueType> probabilitiesToSplitter;
    std::vector<uint64_t> splitterPredecessors;  // states with a non-zero probability to a splitter
    storm::bisimulation::Partition::NonSuperBlockSet nonSuperBlockSet;

    SplitterRefinementCache(storm::bisimulation::Partition const& partition)
        : probabilitiesToSplitter(partition.getNumberOfElements(), storm::utility::zero<ValueType>()), nonSuperBlockSet(partition) {}

    void addProbabilityToSplitter(uint64_t state, ValueType const& probability) {
        STORM_LOG_ASSERT(!storm::utility::isZero(probability), "The probability to add to the splitter must not be zero.");
        if (auto& p = probabilitiesToSplitter[state]; p == storm::utility::zero<ValueType>()) {
            splitterPredecessors.push_back(state);
            p = probability;
        } else {
            p += probability;
        }
    }

    void clear() {
        for (auto const& state : splitterPredecessors) {
            probabilitiesToSplitter[state] = storm::utility::zero<ValueType>();
        }
        splitterPredecessors.clear();
        while (!nonSuperBlockSet.empty()) {
            nonSuperBlockSet.pop();
        }
    }
};

template<typename ValueType>
void refinePartitionBasedOnSplitter(RefinementContext<ValueType>& context, storm::bisimulation::Partition::Block const splitterBlock,
                                    SplitterRefinementCache<ValueType>& cache) {
    auto& partition = context.partition;
    for (auto currentState : splitterBlock) {
        // Compute probability to enter splitter block for each predecessor
        for (const auto& predecessorEntry : context.backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = partition.getBlockOfElement(predecessorState);
            if (predecessorBlock.size() > 1) {  // No need to try to split singleton blocks
                cache.addProbabilityToSplitter(predecessorState, predecessorEntry.getValue());
                cache.nonSuperBlockSet.insert(predecessorBlock);
            }
        }
    }

    while (!cache.nonSuperBlockSet.empty()) {
        auto predecessorBlockToSplit = cache.nonSuperBlockSet.pop();
        // First split the block by whether it is a predecessor of the splitter block or not
        // We do this by either iterating over the splitterPredecessors or the predecessorBlockToSplit, depending on what is shorter.
        auto [noPredecessors, predecessors] =
            cache.splitterPredecessors.size() < predecessorBlockToSplit.size()
                ? partition.splitBlockByRange(predecessorBlockToSplit, cache.splitterPredecessors)
                : partition.splitBlockByPredicate(
                      predecessorBlockToSplit, [&cache](auto const& state) { return !storm::utility::isZero(cache.probabilitiesToSplitter[state]); });

        STORM_LOG_ASSERT(!predecessors.empty(), "The predecessor block should contain at least one predecessor state.");
        bool wasSplit = noPredecessors.size() > 0;

        if (wasSplit) {
            // add the block of states with no transition to the current splitter
            context.queue.insert(noPredecessors);
        }

        if constexpr (storm::IsIntervalType<ValueType>) {
            // TODO: there is no sound tolerance-based order for interval-valued splitter probabilities yet, so we skip the
            // magnitude-based split for interval models for now (predecessors are still split by presence/absence above).
        } else {
            auto const less = [&cache](uint64_t const state1, uint64_t const state2) {
                return cache.probabilitiesToSplitter[state1] < cache.probabilitiesToSplitter[state2];
            };
            if (storm::utility::isZero(context.tolerance)) {
                // Attention: Do not short circuit, i.e., wasSplit = wasSplit || foo() might not execute foo()
                wasSplit |= partition.splitBlockByOrder(predecessors, less);
            } else {
                auto const lessTolerance = [&cache, &context](uint64_t const state1, uint64_t const state2) {
                    return cache.probabilitiesToSplitter[state1] + context.tolerance < cache.probabilitiesToSplitter[state2];
                };
                // Attention: Do not short circuit
                wasSplit |= partition.splitBlockByOrder(predecessors, less, lessTolerance);
            }
        }

        if (wasSplit) {
            // Add all remaining blocks that were split to splitter queue
            context.queue.erase(predecessorBlockToSplit);
            partition.forEachSubBlock(predecessors, [&context](auto const& block) { context.queue.insert(block); });
        }
    }

    // Reset the touched entries of the probabilitiesToCurrentSplitter vector
    cache.clear();
}

}  // namespace detail

template<typename ValueType>
void performPartitionRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                std::optional<std::vector<uint64_t>> const& choiceClasses) {
    detail::RefinementContext<ValueType> context(model, partition, storm::utility::convertNumber<storm::IntervalBaseType<ValueType>>(1e-6));  // todo
    STORM_LOG_ASSERT(model.isNondeterministicModel() || !choiceClasses.has_value(), "Choice classes should only be given for nondeterministic models.");

    // Initially, add all current blocks to the queue.
    partition.forEachBlock([&context](auto const& block) { context.queue.insert(block); });

    // Perform signature-based or splitter-based refinement based on input
    if (model.isNondeterministicModel()) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Signature-based refinement is not supported for nondeterministic models.");
        // idea for signature-based refinement:
        // flag blocks that are known to be stable w.r.t. the current partition
        // Invariant: all unflaged blocks are in the queue (there might also be flagged ones in the queue)
        // take the smallest block B from the queue
        // Split B  into B=B_1 cup B_2 cup ... cup B_n using signature (if B is flagged, we already know that B=B_n=B_1)
        // If refinement of B had no effect, i.e. it  was already stable (B_n=B_1=B) then flag B
        // Compute for each predecessor state s of B an uint64_t s_i with the following interpretation:
        // s_i == 0: s has no transition to B
        // 0 < s_i <= n: s_i has a transition to B_i and no other sub-block of B
        // s_i > n: s_i is a hash encoding the subset of {B_1,...,B_n} that s has transitions to (can enforce this by always setting the most significant bit to
        // true) Note: if B is flagged, then s_i is either 0 or 1 Then, for each predecessor block B' of B:
        //   if B' and B are flagged: continue (assert that  s_i=1 for all s in B')
        //   split B' according to s_i.
        //   if B' is flagged (and B is not flagged):
        //      - there shouldn't be any s in B' with s_i = 0 (assert that! could also make the split of B' above more efficient)
        //      - for subblocks of B' with 0 < s_i <= n, all transitions that lead to B also lead to B_i so the split of B does not introduce any distinguishing
        //           behavior within such blocks. flag those subblocks of B' as well
        //      - subblocks of B' with s_i > n should no longer be flagged (including the case where B' was not split)
        //    if B' was split above: add all subblocks of B' to the queue
        //    else if B' is unflagged: if it was already unflagged, assert that it's already in the queue. Add B' to the queue
        //   Nothing needs to be done if B' is flagged and not split. It might or might not be already in the queue

    } else {
        detail::SplitterRefinementCache<ValueType> refinementCache(partition);

        // Then perform the actual splitting until there are no more splitters.
        while (!context.queue.empty()) {
            auto splitterBlock = *context.queue.begin();
            context.queue.erase(context.queue.begin());
            STORM_LOG_ASSERT(!partition.isProperSuperBlock(splitterBlock), "Broken invariant: the queue should not contain blocks that have been split.");
            refinePartitionBasedOnSplitter(context, splitterBlock, refinementCache);
        }
    }
}

template void performPartitionRefinement<double>(storm::models::sparse::Model<double> const& model, storm::bisimulation::Partition& partition,
                                                 std::optional<std::vector<uint64_t>> const& choiceClasses);
template void performPartitionRefinement<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                storm::bisimulation::Partition& partition,
                                                                std::optional<std::vector<uint64_t>> const& choiceClasses);
template void performPartitionRefinement<storm::RationalFunction>(storm::models::sparse::Model<storm::RationalFunction> const& model,
                                                                  storm::bisimulation::Partition& partition,
                                                                  std::optional<std::vector<uint64_t>> const& choiceClasses);
template void performPartitionRefinement<storm::Interval>(storm::models::sparse::Model<storm::Interval> const& model, storm::bisimulation::Partition& partition,
                                                          std::optional<std::vector<uint64_t>> const& choiceClasses);
template void performPartitionRefinement<storm::RationalInterval>(storm::models::sparse::Model<storm::RationalInterval> const& model,
                                                                  storm::bisimulation::Partition& partition,
                                                                  std::optional<std::vector<uint64_t>> const& choiceClasses);

}  // namespace storm::bisimulation

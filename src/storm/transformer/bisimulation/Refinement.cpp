#include "storm/transformer/bisimulation/Refinement.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
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
        : model(model), partition(partition), backwardTransitions(model.getBackwardTransitions()), comparator(tolerance) {}

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition& partition;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    storm::utility::ConstantsComparator<storm::IntervalBaseType<ValueType>> const comparator;
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
            // todo: use a second order in this call: one for initial sort, one for split
            // Attention: Do not short circuit, i.e., wasSplit = wasSplit || foo() might not execute foo()
            wasSplit |= partition.splitBlockByOrder(predecessors, [&cache, &context](auto const& a, auto const& b) {
                return context.comparator.isLess(cache.probabilitiesToSplitter[a], cache.probabilitiesToSplitter[b]);
            });
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
void performPartitionRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition) {
    detail::RefinementContext<ValueType> context(model, partition, storm::utility::convertNumber<storm::IntervalBaseType<ValueType>>(1e-6));  // todo

    // Initially, add all current blocks to the queue.
    partition.forEachBlock([&context](auto const& block) { context.queue.insert(block); });

    detail::SplitterRefinementCache<ValueType> refinementCache(partition);

    // Then perform the actual splitting until there are no more splitters.
    while (!context.queue.empty()) {
        auto splitterBlock = *context.queue.begin();
        context.queue.erase(context.queue.begin());
        STORM_LOG_ASSERT(!partition.isProperSuperBlock(splitterBlock), "Broken invariant: the queue should not contain blocks that have been split.");
        refinePartitionBasedOnSplitter(context, splitterBlock, refinementCache);
    }
}

template void performPartitionRefinement<double>(storm::models::sparse::Model<double> const& model, storm::bisimulation::Partition& partition);
template void performPartitionRefinement<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                storm::bisimulation::Partition& partition);
template void performPartitionRefinement<storm::RationalFunction>(storm::models::sparse::Model<storm::RationalFunction> const& model,
                                                                  storm::bisimulation::Partition& partition);
template void performPartitionRefinement<storm::Interval>(storm::models::sparse::Model<storm::Interval> const& model,
                                                          storm::bisimulation::Partition& partition);
template void performPartitionRefinement<storm::RationalInterval>(storm::models::sparse::Model<storm::RationalInterval> const& model,
                                                                  storm::bisimulation::Partition& partition);

}  // namespace storm::bisimulation

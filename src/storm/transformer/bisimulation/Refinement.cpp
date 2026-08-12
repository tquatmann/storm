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
    RefinementContext(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                      storm::bisimulation::Partition& partition, storm::IntervalBaseType<ValueType> tolerance)
        : model(model), choiceClasses(choiceClasses), partition(partition), backwardTransitions(model.getBackwardTransitions()), tolerance(tolerance) {}

    storm::models::sparse::Model<ValueType> const& model;
    std::optional<std::vector<uint64_t>> const& choiceClasses;
    storm::bisimulation::Partition& partition;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    ValueType const tolerance;
    storm::bisimulation::Partition::OrderedBlockSet queue;
};

template<typename ValueType>
struct SplitterRefinementCache {
    std::vector<ValueType> probabilitiesToSplitter;
    std::vector<uint64_t> splitterPredecessorStates;  // states with a non-zero probability to a splitter
    Partition::NonSuperBlockSet splitterPredecessorBlocks;

    explicit SplitterRefinementCache(Partition const& partition)
        : probabilitiesToSplitter(partition.getNumberOfElements(), storm::utility::zero<ValueType>()), splitterPredecessorBlocks(partition) {}

    void addProbabilityToSplitter(uint64_t state, ValueType const& probability) {
        STORM_LOG_ASSERT(!storm::utility::isZero(probability), "The probability to add to the splitter must not be zero.");
        if (auto& p = probabilitiesToSplitter[state]; p == storm::utility::zero<ValueType>()) {
            splitterPredecessorStates.push_back(state);
            p = probability;
        } else {
            p += probability;
        }
    }

    void clear() {
        for (auto const& state : splitterPredecessorStates) {
            probabilitiesToSplitter[state] = storm::utility::zero<ValueType>();
        }
        splitterPredecessorStates.clear();
        while (!splitterPredecessorBlocks.empty()) {
            splitterPredecessorBlocks.pop();
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
                cache.splitterPredecessorBlocks.insert(predecessorBlock);
            }
        }
    }

    while (!cache.splitterPredecessorBlocks.empty()) {
        auto predecessorBlockToSplit = cache.splitterPredecessorBlocks.pop();
        // First split the block by whether a state is a predecessor of the splitter block or not
        // We do this by either iterating over the splitterPredecessors or the predecessorBlockToSplit, depending on what is shorter.
        auto [noPredecessors, predecessors] =
            cache.splitterPredecessorStates.size() < predecessorBlockToSplit.size()
                ? partition.splitBlockByRange(predecessorBlockToSplit, cache.splitterPredecessorStates)
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

template<typename ValueType>
struct ChoiceSignature {
    uint64_t const choiceClass;
    Partition::OrderedBlockMap<ValueType> distr;  // todo: try alternatives

    ChoiceSignature(uint64_t const choiceIndex, RefinementContext<ValueType> const& context)
        : choiceClass(context.choiceClasses ? (*context.choiceClasses)[choiceIndex] : 0) {
        for (auto const& entry : context.model.getTransitionMatrix().getRow(choiceIndex)) {
            if (!storm::utility::isZero(entry.getValue())) {
                auto emplace_res = distr.emplace(context.partition.getBlockOfElement(entry.getColumn()), entry.getValue());
                if (!emplace_res.second) {
                    emplace_res.first->second += entry.getValue();
                }
            }
        }
    }
    // TODO: handle tolerance!
    bool operator<(ChoiceSignature const& other) const {
        if (choiceClass != other.choiceClass) {
            return choiceClass < other.choiceClass;
        }
        // std::map's relational operators are unavailable here since Block (a std::span) has no operator</operator==,
        // so we compare manually using BlockCompare (both maps are already ordered by it).
        Partition::BlockCompare const blockLess;
        auto it1 = distr.begin();
        auto it2 = other.distr.begin();
        for (; it1 != distr.end() && it2 != other.distr.end(); ++it1, ++it2) {
            if (blockLess(it1->first, it2->first)) {
                return true;
            }
            if (blockLess(it2->first, it1->first)) {
                return false;
            }
            if (it1->second != it2->second) {
                return it1->second < it2->second;
            }
        }
        return it1 == distr.end() && it2 != other.distr.end();
    }

    bool operator==(ChoiceSignature const& other) const {
        if (choiceClass != other.choiceClass || distr.size() != other.distr.size()) {
            return false;
        }
        Partition::BlockCompare const blockLess;
        auto it2 = other.distr.begin();
        for (auto it1 = distr.begin(); it1 != distr.end(); ++it1, ++it2) {
            if (blockLess(it1->first, it2->first) || blockLess(it2->first, it1->first) || it1->second != it2->second) {
                return false;
            }
        }
        return true;
    };
};

template<typename ValueType>
struct StateSignature {
    std::set<ChoiceSignature<ValueType>> choices;
    StateSignature() = default;
    StateSignature(uint64_t const stateIndex, RefinementContext<ValueType> const& context) {
        for (uint64_t const choiceIndex : context.model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
            choices.emplace(choiceIndex, context);
        }
    }

    bool operator<(StateSignature const& other) const {
        return choices < other.choices;
    }
};

template<typename ValueType>
struct SignatureRefinementCache {
    std::vector<StateSignature<ValueType>> stateSignatures;                 // the signature of each state
    storm::bisimulation::Partition::OrderedBlockSet stableCandidateBlocks;  // Blocks that might be considered stable w.r.t. the current partition
    std::vector<std::set<uint64_t>> reachablePivotBlocks;                   // used to store which sub-blocks of the pivot can be reached for each state
    std::vector<uint64_t> pivotPredecessorStates;                           // states with a non-zero probability to the pivot block
    Partition::NonSuperBlockSet pivotPredecessorBlocks;                     // blocks with a non-zero probability to the pivot block

    explicit SignatureRefinementCache(storm::bisimulation::Partition const& partition)
        : stateSignatures(partition.getNumberOfElements()),
          stableCandidateBlocks(),
          reachablePivotBlocks(partition.getNumberOfElements()),
          pivotPredecessorBlocks(partition) {}

    void addReachablePivotBlock(uint64_t const state, uint64_t const localSubBlockIndex) {
        auto& blocks = reachablePivotBlocks[state];
        if (blocks.empty()) {
            pivotPredecessorStates.push_back(state);
        }
        blocks.insert(localSubBlockIndex);
    }

    void clearPivotPredecessorData() {
        for (uint64_t const state : pivotPredecessorStates) {
            reachablePivotBlocks[state].clear();
        }
        pivotPredecessorStates.clear();
        while (!pivotPredecessorBlocks.empty()) {
            pivotPredecessorBlocks.pop();
        }
    }
};

template<typename ValueType>
void refinePartitionBasedOnSignature(RefinementContext<ValueType>& context, storm::bisimulation::Partition::Block const pivotBlock,
                                     SignatureRefinementCache<ValueType>& cache) {
    // Split the pivot block B into B=B_1 cup B_2 cup ... cup B_n using signature refinement
    // First update the state signatures
    for (uint64_t const state : pivotBlock) {
        cache.stateSignatures[state] = StateSignature<ValueType>(state, context);
    }
    // Then perform the signature-based split
    bool const pivotStable = !context.partition.splitBlockByOrder(
        pivotBlock, [&cache](uint64_t const state1, uint64_t const state2) { return cache.stateSignatures[state1] < cache.stateSignatures[state2]; });

    if (pivotStable && cache.stableCandidateBlocks.contains(pivotBlock)) {
        // When the current pivot block is stable, we can continue with the next pivot
        cache.stableCandidateBlocks.erase(pivotBlock);
        return;
    }
    // When the pivot block has been split, all the predecessor blocks need to checked again.
    // While we could just add all those predecessors to the queue, we try to split them first based on simple, graph-based criteria so that (expensive)
    // signature refinement is hopefully only applied to smaller blocks. Specifically, we split predecessor blocks based on which set of sub-blocks of the pivot
    // they can reach.

    // Gather predecessors of the pivot and their reachable sub-blocks
    uint64_t subBlockIndex = 0;
    context.partition.forEachSubBlock(pivotBlock, [&context, &cache, &subBlockIndex](auto const& subBlock) {
        for (uint64_t const state : subBlock) {
            for (auto const& predecessorEntry : context.backwardTransitions.getRow(state)) {
                auto predecessorState = predecessorEntry.getColumn();
                auto predecessorBlock = context.partition.getBlockOfElement(predecessorState);
                if (predecessorBlock.size() > 1) {
                    // No need to investigate singleton predecessor blocks as they cannot be split any further.
                    cache.addReachablePivotBlock(predecessorState, subBlockIndex);
                    cache.pivotPredecessorBlocks.insert(predecessorBlock);
                }
            }
        }
        ++subBlockIndex;
    });
    // Apply splitting of the pivot predecessors (similar to splitter-based refinement)
    while (!cache.pivotPredecessorBlocks.empty()) {
        auto const predecessorBlock = cache.pivotPredecessorBlocks.pop();
        // Split the predecessor block according to which sub-blocks of the pivot-block can be reached.

        // Split the block by whether a state is a predecessor of the pivot block or not
        // We do this by either iterating over the pivotPredecessorStates or the predecessorBlock, depending on what is shorter.
        auto [noPredecessors, predecessors] =
            cache.pivotPredecessorStates.size() < predecessorBlock.size()
                ? context.partition.splitBlockByRange(predecessorBlock, cache.pivotPredecessorStates)
                : context.partition.splitBlockByPredicate(predecessorBlock, [&cache](auto const& state) { return !cache.reachablePivotBlocks[state].empty(); });

        // At least one state should be a predecessor of the pivot block (otherwise we wouldn't have found that block above)
        STORM_LOG_ASSERT(!predecessors.empty(), "The predecessor block should contain at least one predecessor state.");

        // Now apply the splitting based on which sub-blocks of the pivot block can be reached.
        // Attention: Do not short circuit, i.e., wasSplit = wasSplit || foo() might not execute foo()
        context.partition.splitBlockByOrder(predecessors, [&cache](uint64_t const state1, uint64_t const state2) {
            return cache.reachablePivotBlocks[state1] < cache.reachablePivotBlocks[state2];
        });

        if (context.partition.isProperSuperBlock(predecessorBlock)) {
            // Erase the super block from the queue and the stable candidates (it might or might not be in there)
            context.queue.erase(predecessorBlock);
            cache.stableCandidateBlocks.erase(predecessorBlock);
            // Add all sub-blocks to the queue
            context.partition.forEachSubBlock(predecessorBlock, [&context](auto const& block) { context.queue.insert(block); });
        } else {
            // The simple, graph-based splitting was not effective. We must add the entire predecessorBlock to the queue
            cache.stableCandidateBlocks.insert(predecessorBlock);
            context.queue.insert(predecessorBlock);
        }
    }
}

}  // namespace detail

template<typename ValueType>
void performPartitionRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                std::optional<std::vector<uint64_t>> const& choiceClasses) {
    detail::RefinementContext<ValueType> context(model, choiceClasses, partition,
                                                 storm::utility::convertNumber<storm::IntervalBaseType<ValueType>>(1e-6));  // todo
    STORM_LOG_ASSERT(model.isNondeterministicModel() || !choiceClasses.has_value(), "Choice classes should only be given for nondeterministic models.");

    // Initially, add all current blocks to the queue.
    partition.forEachBlock([&context](auto const& block) { context.queue.insert(block); });

    // Perform signature-based or splitter-based refinement based on input
    if (model.isNondeterministicModel()) {
        detail::SignatureRefinementCache<ValueType> cache(partition);
        while (!context.queue.empty()) {
            // take the smallest block from the queue
            auto const pivotBlock = *context.queue.begin();
            context.queue.erase(context.queue.begin());
            STORM_LOG_ASSERT(!partition.isProperSuperBlock(pivotBlock), "Broken invariant: the queue should not contain blocks that have been split.");
            // Split the pivotBlock based on its signature and split the predecessor blocks based on a simple, graph-based criterion
            detail::refinePartitionBasedOnSignature(context, pivotBlock, cache);
        }
    } else {
        detail::SplitterRefinementCache<ValueType> refinementCache(partition);

        // Perform the splitting until there are no more splitters.
        while (!context.queue.empty()) {
            // take the smallest block from the queue
            auto const splitterBlock = *context.queue.begin();
            context.queue.erase(context.queue.begin());
            STORM_LOG_ASSERT(!partition.isProperSuperBlock(splitterBlock), "Broken invariant: the queue should not contain blocks that have been split.");
            // Split the predecessor blocks and add them to the queue
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

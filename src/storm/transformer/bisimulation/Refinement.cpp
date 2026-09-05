#include "storm/transformer/bisimulation/Refinement.h"

#include <ranges>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Signatures.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

namespace detail {

/*!
 * Stores a mapping from states to ValueType
 * Adding values, reading values, and clearing can all be done in constant time.
 * However, the size of the map is linear in the total number of states
 */
template<typename ValueType>
class StateMapping {
   public:
    explicit StateMapping(uint64_t const numStates) : values(numStates, defaultValue()) {}

    /*!
     * @return Retrieves the currently stored values
     */
    std::vector<ValueType> const& getValues() const {
        return values;
    }

    /*!
     * @return the list of states currently holding a non-zero value
     */
    std::vector<uint64_t> const& getNonDefaultStates() const {
        return nonDefaultStates;
    }

    /*!
     * Adds value to the currently mapped value of the given state
     */
    template<typename T>
    void addValue(uint64_t const state, T value) {
        if constexpr (std::is_same_v<ValueType, std::set<uint64_t>>) {
            if (values[state].empty()) {
                nonDefaultStates.push_back(state);
            }
            values[state].insert(value);
        } else {
            STORM_LOG_ASSERT(!storm::utility::isZero(value), "Did not expect adding 0 probability");
            if (storm::utility::isZero(values[state])) {
                nonDefaultStates.push_back(state);
                values[state] = std::move(value);
            } else {
                values[state] += std::move(value);
            }
        }
    }

    /*!
     * Clears the set, i.e., writes the default value for all states.
     */
    void clear() {
        for (auto const& state : nonDefaultStates) {
            values[state] = defaultValue();
        }
        nonDefaultStates.clear();
    }

   private:
    static ValueType defaultValue() {
        if constexpr (std::is_same_v<ValueType, std::set<uint64_t>>) {
            return {};  // empty set
        } else {
            return storm::utility::zero<ValueType>();
        }
    }

    std::vector<ValueType> values;           // stores the value for each state
    std::vector<uint64_t> nonDefaultStates;  // stores those states with a non-default value
};

template<typename ValueType>
struct SplitterRefinementContext {
    SplitterRefinementContext(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition, ValueType const tolerance)
        : model(model), partition(partition), backwardTransitions(model.getBackwardTransitions()), tolerance(tolerance), cache(partition) {}

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition& partition;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    ValueType const tolerance;
    storm::bisimulation::Partition::OrderedBlockSet queue;

    struct Cache {
        Cache(Partition const& partition) : predecessorToSplitterProbabilities(partition.getNumberOfElements()), predecessorBlocks(partition) {}
        StateMapping<ValueType> predecessorToSplitterProbabilities;
        Partition::NonSuperBlockSet predecessorBlocks;
    } cache;
};

template<typename ValueType>
void refinePartitionBasedOnSplitter(SplitterRefinementContext<ValueType>& context, storm::bisimulation::Partition::Block const splitterBlock) {
    auto& predecessorToSplitterProbabilities = context.cache.predecessorToSplitterProbabilities;
    auto& predecessorBlocks = context.cache.predecessorBlocks;

    for (auto currentState : splitterBlock) {
        // Compute probability to enter splitter block for each predecessor
        for (const auto& predecessorEntry : context.backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = context.partition.getBlockOfElement(predecessorState);
            if (predecessorBlock.size() > 1) {  // No need to try to split singleton blocks
                predecessorToSplitterProbabilities.addValue(predecessorState, predecessorEntry.getValue());
                predecessorBlocks.insert(predecessorBlock);
            }
        }
    }

    while (!predecessorBlocks.empty()) {
        auto predecessorBlockToSplit = predecessorBlocks.pop();
        // First split the block by whether a state is a predecessor of the splitter block or not
        // We do this by either iterating over the splitterPredecessors or the predecessorBlockToSplit, depending on what is shorter.
        auto [noPredecessors, predecessors] =
            predecessorToSplitterProbabilities.getNonDefaultStates().size() < predecessorBlockToSplit.size()
                ? context.partition.splitBlockByRange(predecessorBlockToSplit, predecessorToSplitterProbabilities.getNonDefaultStates())
                : context.partition.splitBlockByPredicate(
                      predecessorBlockToSplit,
                      [&predecessorToSplitterProbabilities](auto const& state) {
                          return !storm::utility::isZero(predecessorToSplitterProbabilities.getValues()[state]);
                      });

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
            auto const& toSplitterProbs = predecessorToSplitterProbabilities.getValues();
            auto const less = [&toSplitterProbs](uint64_t const state1, uint64_t const state2) { return toSplitterProbs[state1] < toSplitterProbs[state2]; };
            if (storm::utility::isZero(context.tolerance)) {
                // Attention: Do not short circuit, i.e., wasSplit = wasSplit || foo() might not execute foo()
                wasSplit |= context.partition.splitBlockByOrder(predecessors, less);
            } else {
                auto const lessTolerance = [&toSplitterProbs, &context](uint64_t const state1, uint64_t const state2) {
                    return toSplitterProbs[state1] + context.tolerance < toSplitterProbs[state2];
                };
                // Attention: Do not short circuit
                wasSplit |= context.partition.splitBlockByOrder(predecessors, less, lessTolerance);
            }
        }

        if (wasSplit) {
            // Add all remaining blocks that were split to splitter queue
            context.queue.erase(predecessorBlockToSplit);
            context.partition.forEachSubBlock(predecessors, [&context](auto const& block) { context.queue.insert(block); });
        }
    }

    // Reset the predecessorToSplitterProbabilities for the next iteration.
    predecessorToSplitterProbabilities.clear();
}

template<typename ValueType, SignatureMode SignatureMode>
struct SignatureRefinementContext {
    SignatureRefinementContext(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                               Signatures<ValueType, SignatureMode>& signatures)
        : model(model), partition(partition), signatures(signatures), backwardTransitions(model.getBackwardTransitions()), cache(partition) {}

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition& partition;
    storm::bisimulation::Signatures<ValueType, SignatureMode>& signatures;
    storm::storage::SparseMatrix<ValueType> const backwardTransitions;
    storm::bisimulation::Partition::OrderedBlockMap<bool>
        queue;  // stores an extra flag for each element in the queue. The flag indicates whether we enforce exploring the predecessors of the block

    struct Cache {
        Cache(Partition const& partition) : predecessorToPivotBlocks(partition.getNumberOfElements()), predecessorBlocks(partition) {}
        StateMapping<std::set<uint64_t>> predecessorToPivotBlocks;
        Partition::NonSuperBlockSet predecessorBlocks;
    } cache;
};

template<typename ValueType, SignatureMode SignatureMode>
void refinePartitionBasedOnSignature(SignatureRefinementContext<ValueType, SignatureMode>& context, storm::bisimulation::Partition::Block const pivotBlock,
                                     bool const enforcePredecessorExploration) {
    // Split the pivot block B into B=B_1 cup B_2 cup ... cup B_n using signature refinement
    // First update the state signatures
    for (uint64_t const state : pivotBlock) {
        context.signatures.updateStateSignature(state);
    }
    // Then perform the signature-based split
    bool pivotHasBeenSplit = context.partition.splitBlockByOrder(pivotBlock, context.signatures.getSplitOrder(), context.signatures.getSplitCondition());

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
                : context.partition.splitBlockByPredicate(predecessorBlock, [&toPivotBlocks](auto const& state) { return !toPivotBlocks[state].empty(); });

        // At least one state should be a predecessor of the pivot block (otherwise we wouldn't have found that block above)
        STORM_LOG_ASSERT(!predecessors.empty(), "The predecessor block should contain at least one predecessor state.");

        // Now apply the splitting based on which sub-blocks of the pivot block can be reached.
        // If we did not actually split the pivot block, this operation would have no effect.
        if (pivotHasBeenSplit) {
            context.partition.splitBlockByOrder(
                predecessors, [&toPivotBlocks](uint64_t const state1, uint64_t const state2) { return toPivotBlocks[state1] < toPivotBlocks[state2]; });
        } else {
            STORM_LOG_ASSERT(
                std::all_of(predecessors.begin(), predecessors.end(),
                            [&toPivotBlocks](uint64_t const& state) { return toPivotBlocks[state].size() == 1 && *toPivotBlocks[state].begin() == 0; }),
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

template<typename ValueType>
void performSplitterBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                    ValueType const tolerance) {
    static_assert(!storm::IsIntervalType<ValueType>, "Interval types are not supported for splitter-based refinement.");
    // Refinement for interval models requires limiting signatures to feasible intervals. This is rather difficult in a splitter-based setting.
    STORM_LOG_THROW(!model.isNondeterministicModel(), storm::exceptions::InvalidArgumentException,
                    "Splitter-based refinement is only supported for deterministic models.");
    STORM_LOG_THROW((storm::utility::isZero(tolerance) || !std::is_same_v<ValueType, storm::RationalFunction>), storm::exceptions::InvalidArgumentException,
                    "Splitter-based refinement with non-zero tolerance does not apply to parametric models.");
    detail::SplitterRefinementContext<ValueType> context(model, partition, tolerance);
    // Initially, add all current blocks to the queue.
    partition.forEachBlock([&context](auto const& block) { context.queue.insert(block); });

    // Perform the splitting until there are no more splitters.
    while (!context.queue.empty()) {
        // take the smallest block from the queue
        auto const splitterBlock = *context.queue.begin();
        context.queue.erase(context.queue.begin());
        STORM_LOG_ASSERT(!partition.isProperSuperBlock(splitterBlock), "Broken invariant: the queue should not contain blocks that have been split.");
        // Split the predecessor blocks and add them to the queue
        detail::refinePartitionBasedOnSplitter(context, splitterBlock);
    }
}

template<typename ValueType, SignatureMode SignatureMode>
void performSignatureBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                     Signatures<ValueType, SignatureMode>& signatures) {
    static_assert(!storm::IsIntervalType<ValueType>, "Interval types are not yet supported for signature-based refinement.");
    detail::SignatureRefinementContext<ValueType, SignatureMode> context(model, partition, signatures);
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
        }
    });
}

template void performSplitterBasedRefinement<double>(storm::models::sparse::Model<double> const& model, storm::bisimulation::Partition& partition,
                                                     double const tolerance);
template void performSplitterBasedRefinement<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                    storm::bisimulation::Partition& partition, storm::RationalNumber const tolerance);
template void performSplitterBasedRefinement<storm::RationalFunction>(storm::models::sparse::Model<storm::RationalFunction> const& model,
                                                                      storm::bisimulation::Partition& partition, storm::RationalFunction const tolerance);

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

}  // namespace storm::bisimulation

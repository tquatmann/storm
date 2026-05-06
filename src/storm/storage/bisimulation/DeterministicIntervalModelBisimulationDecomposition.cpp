#include "DeterministicIntervalModelBisimulationDecomposition.h"
#include <algorithm>
#include <boost/iterator/zip_iterator.hpp>
#include <chrono>
#include <iomanip>
#include <storm/models/sparse/Dtmc.h>
#include <unordered_map>
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"

#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/graph.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::DeterministicIntervalModelBisimulationDecomposition(
    const ModelType& model, const typename BisimulationDecomposition<ModelType>::Options& options)
    : storm::storage::BisimulationDecomposition<ModelType>(model, options),
      probabilitiesToCurrentSplitter(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      probabilitiesToOtherBlocks(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      touchedProbabilitiesToSplitter(model.getNumberOfStates(), false),
      originalStateDistributions(model.getNumberOfStates()) {}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::initializeMeasureDrivenPartition() {
    // TODO: implement actual measure-driven partitioning
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::initializeLabelBasedPartition() {
    BisimulationDecomposition<ModelType>::initializeLabelBasedPartition();
    postProcessInitialPartition();
}

template<typename ModelType>
Distribution<typename ModelType::ValueType, sparse::state_type> DeterministicIntervalModelBisimulationDecomposition<ModelType>::getClampedDistribution(
    Distribution<ValueType, storm::storage::sparse::state_type> distribution) const {
    auto clampedDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
    clampedDistribution.reserve(distribution.size());

    for (auto it = distribution.begin(); it != distribution.end(); ++it) {
        auto key = it->first;
        auto interval = it->second;

        clampedDistribution.addProbability(key, clampIntervalToProbabilisticInterval(interval));
    }

    return clampedDistribution;
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(
    storm::storage::Partition::Block splitterBlock, std::deque<typename bisimulation::Partition::Block>& splitterQueue,
    bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::bisimulation::Partition::BlockSet blocksToSplit;

    for (auto currentState : splitterBlock) {
        for (const auto& predecessorEntry : this->backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);
            auto intervalTransitionProbability = predecessorEntry.getValue();

            if (!possiblyNeedsRefinement(predecessorBlock)) {
                continue;
            }

            // Compute interval going to the splitter, i.e., I(s, C), and initialize interval leading to all other blocks, i.e., I(s, \Pi \setminus C)
            if (touchedProbabilitiesToSplitter.get(predecessorState)) {
                probabilitiesToCurrentSplitter[predecessorState] += intervalTransitionProbability;
            } else {
                probabilitiesToCurrentSplitter[predecessorState] = intervalTransitionProbability;
                probabilitiesToOtherBlocks[predecessorState] = storm::utility::zero<ValueType>();
                touchedProbabilitiesToSplitter.set(predecessorState, true);
            }

            // Remember which blocks contain predecessors to split w.r.t. the splitter afterwards
            blocksToSplit.insert(predecessorBlock);
        }
    }

    // Compute interval to all blocks except the splitter, i.e., I(s, \Pi \setminus C)
    for (auto predecessorState : touchedProbabilitiesToSplitter) {
        for (const auto& entry : this->model.getTransitionMatrix().getRow(predecessorState)) {
            auto targetState = entry.getColumn();

            if (!this->partition.contains(splitterBlock, targetState)) {
                probabilitiesToOtherBlocks[predecessorState] += entry.getValue();
            }
        }
    }

    for (auto predecessorBlockToSplit : blocksToSplit) {
        // First split the block by whether it is a predecessor of the splitter block or not
        auto [noPredecessors, predecessors] =
            this->partition.splitBlockByPredicate(predecessorBlockToSplit, [this](auto const& state) { return touchedProbabilitiesToSplitter.get(state); });

        if (noPredecessors.size() > 0) {
            if (!enqueuedSplitterBlocks.contains(noPredecessors)) {
                splitterQueue.push_back(noPredecessors);
                enqueuedSplitterBlocks.insert(noPredecessors);
            }
        }

        if (predecessors.size() > 0) {
            bool wasSplit = this->partition.splitBlockByOrder(predecessors, [this](auto const& a, auto const& b) {
                if (!touchedProbabilitiesToSplitter.get(a) || !touchedProbabilitiesToSplitter.get(b)) {
                    STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Comparing states that were not touched: " << a << ", " << b);
                }

                auto feasibleIntervalOfStateA =
                    computeFeasibleIntervalBasedOnAggregatedIntervals(probabilitiesToCurrentSplitter[a], probabilitiesToOtherBlocks[a]);
                auto feasibleIntervalOfStateB =
                    computeFeasibleIntervalBasedOnAggregatedIntervals(probabilitiesToCurrentSplitter[b], probabilitiesToOtherBlocks[b]);

                // Compare the "raw" interval bounds
                // return probabilitiesToCurrentSplitter[a].lower() < probabilitiesToCurrentSplitter[b].lower() ||
                //        (probabilitiesToCurrentSplitter[a].lower() == probabilitiesToCurrentSplitter[b].lower() &&
                //         probabilitiesToCurrentSplitter[a].upper() < probabilitiesToCurrentSplitter[b].upper());

                // Compare the feasible intervals
                return feasibleIntervalOfStateA.lower() < feasibleIntervalOfStateB.lower() ||
                       (feasibleIntervalOfStateA.lower() == feasibleIntervalOfStateB.lower() &&
                        feasibleIntervalOfStateA.upper() < feasibleIntervalOfStateB.upper());
            });

            if (wasSplit) {
                this->partition.forEachSubBlock(predecessors, [&splitterQueue, &enqueuedSplitterBlocks](auto const& block) {
                    if (!enqueuedSplitterBlocks.contains(block)) {
                        splitterQueue.push_back(block);
                        enqueuedSplitterBlocks.insert(block);
                    }
                });
            }
        }
    }

    touchedProbabilitiesToSplitter.clear();
}

template<typename ModelType>
ModelType::ValueType DeterministicIntervalModelBisimulationDecomposition<ModelType>::computeFeasibleIntervalBasedOnAggregatedIntervals(
    ValueType intervalToSplitter, ValueType intervalToOtherBlocks) const {
    // Normalize intervals
    ValueType normalizedIntervalToSplitter = clampIntervalToProbabilisticInterval(intervalToSplitter);
    ValueType normalizedIntervalToOtherBlocks = clampIntervalToProbabilisticInterval(intervalToOtherBlocks);

    // Compute feasible interval
    storm::IntervalBaseType<ValueType> lowerBoundOfFeasibleInterval;
    if (normalizedIntervalToSplitter.lower() > storm::utility::one<IntervalBaseType<ValueType>>() - normalizedIntervalToOtherBlocks.upper()) {
        lowerBoundOfFeasibleInterval = normalizedIntervalToSplitter.lower();
    } else {
        lowerBoundOfFeasibleInterval = storm::utility::one<IntervalBaseType<ValueType>>() - normalizedIntervalToOtherBlocks.upper();
    }

    storm::IntervalBaseType<ValueType> upperBoundOfFeasibleInterval;
    if (normalizedIntervalToSplitter.upper() < storm::utility::one<IntervalBaseType<ValueType>>() - normalizedIntervalToOtherBlocks.lower()) {
        upperBoundOfFeasibleInterval = normalizedIntervalToSplitter.upper();
    } else {
        upperBoundOfFeasibleInterval = storm::utility::one<IntervalBaseType<ValueType>>() - normalizedIntervalToOtherBlocks.lower();
    }

    return ValueType(lowerBoundOfFeasibleInterval, carl::BoundType::WEAK, upperBoundOfFeasibleInterval, carl::BoundType::WEAK);
}

template<typename ModelType>
bool DeterministicIntervalModelBisimulationDecomposition<ModelType>::possiblyNeedsRefinement(std::span<uint64_t const> block) const {
    return block.size() > 1 && !this->absorbingBlocks.contains(block.front());
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::refineBlockBasedOnEpsilonSignature(
    std::span<uint64_t const> block, std::deque<typename bisimulation::Partition::Block>& blocksQueue, bisimulation::Partition::BlockSet& enqueuedBlocks,
    double epsilon) {
    // Cluster states of this block in groups
    std::vector<std::pair<std::vector<storm::storage::sparse::state_type>, storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>>> groups;

    for (auto s : block) {
        bool placed = false;
        for (auto& g : groups) {
            bool fits = true;
            auto candidateDistribution = computeEnhancedDistribution(g.second, originalStateDistributions[s]);

            // Check if enhancement satisfies \varepsilon-criterion via our upper bound w.r.t. s itself
            if (2 * computeDeltaForState(originalStateDistributions[s], candidateDistribution) > epsilon) {
                continue;
            }

            for (auto t : g.first) {
                // Compute delta between the state distribution of t and the candidate distribution based on the current group distribution including the
                // interval extensions based on candidate state t
                auto delta = computeDeltaForState(originalStateDistributions[t], candidateDistribution);
                if (2 * delta > epsilon) {
                    fits = false;
                    break;
                }
            }

            // Add state to current group
            if (fits) {
                g.first.push_back(s);
                g.second = candidateDistribution;
                placed = true;
                break;
            }
        }

        // If the state was not placed yet, create its own group
        if (!placed) {
            groups.push_back(
                std::pair<std::vector<storm::storage::sparse::state_type>, storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>>(
                    {s}, originalStateDistributions[s]));
        }
    }

    if (groups.size() <= 1)
        return;  // Nothing to split

    // Split block by grouping
    bool wasSplit = this->partition.splitBlockByOrder(block, [&](auto a, auto b) {
        int groupOfA = -1, groupOfB = -1;
        for (int i = 0; i < groups.size(); ++i) {
            if (std::find(groups[i].first.begin(), groups[i].first.end(), a) != groups[i].first.end())
                groupOfA = i;
            if (std::find(groups[i].first.begin(), groups[i].first.end(), b) != groups[i].first.end())
                groupOfB = i;
        }
        return (groupOfA < groupOfB);
    });

    if (wasSplit) {
        // Place blocks of predecessors into queue
        this->partition.forEachSubBlock(block, [this, &blocksQueue, &enqueuedBlocks](auto const& subBlock) {
            for (auto state : subBlock) {
                for (auto& transition : this->backwardTransitions.getRow(state)) {
                    auto predecessorState = transition.getColumn();
                    auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);

                    // Place predecessor block into queue only if it is not already there
                    if (predecessorBlock.size() > 1 && !enqueuedBlocks.contains(predecessorBlock)) {
                        blocksQueue.push_back(predecessorBlock);
                        enqueuedBlocks.insert(predecessorBlock);
                    }

                    // Recompute distribution of predecessor state
                    auto distribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
                    for (auto const& e : this->model.getTransitionMatrix().getRow(predecessorState)) {
                        distribution.addProbability(this->partition.getBlockOfElement(e.getColumn()).front(), e.getValue());
                    }
                    originalStateDistributions[predecessorState] = getClampedDistribution(distribution);
                }
            }
        });
    }
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::buildQuotient() {
    // if (this->options.getUsesEpsilon()) {
    //     debugFindMergeableFinalBlocks(this->options.getEpsilon());
    //     validateEpsilonStability(this->options.getEpsilon());
    // } else {
    //     checkCurrentPartitionByExactFeasibleIntervals();
    // }

    // In order to create the quotient model, we need to construct
    // (a) the new transition matrix,
    // (b) the new labeling,
    // (c) the new reward structures.

    // Prepare a matrix builder for (a).
    storm::storage::SparseMatrixBuilder<ValueType> builder(this->partition.getNumberOfBlocks(), this->partition.getNumberOfBlocks());

    // Prepare the new state labeling for (b).
    storm::models::sparse::StateLabeling newLabeling(this->partition.getNumberOfBlocks());
    std::set<std::string> atomicPropositionsSet = this->options.respectedAtomicPropositions.value();
    atomicPropositionsSet.insert("init");
    std::vector<std::string> atomicPropositions = std::vector<std::string>(atomicPropositionsSet.begin(), atomicPropositionsSet.end());
    for (auto const& ap : atomicPropositions) {
        newLabeling.addLabel(ap);
    }

    // If the model had state rewards, we need to build the state rewards for the quotient as well.
    std::optional<std::vector<ValueType>> stateRewards;
    if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
        stateRewards = std::vector<ValueType>(this->partition.getNumberOfBlocks());
    }

    // Now build (a) and (b) by traversing all blocks.

    // TODO: create mapping from representative state to unique identifier
    std::map<uint64_t, uint64_t> blocksMapping;

    auto blockIndex = 0;
    this->partition.forEachBlock([&](auto const& block) {
        blocksMapping.emplace(block.front(), blockIndex);
        blockIndex++;
    });

    this->partition.forEachBlock([&](auto const& block) {
        // Pick one representative state. For strong bisimulation it doesn't matter which state it is, because
        // they all behave equally.
        auto representativeState = block.front();

        // TODO: Handle weak bisimulation case (non-silent)

        blockIndex = blocksMapping.at(representativeState);

        // If the block is absorbing, we simply add a self-loop.
        if (this->absorbingBlocks.contains(representativeState)) {
            builder.addNextValue(blockIndex, blockIndex, storm::utility::one<ValueType>());

            // If the block has a special representative state, we retrieve it now.
            representativeState = this->absorbingBlocks.at(representativeState);

            // Add all selected atomic propositions that hold in the representative state to the state
            // representing the block.
            for (auto const& ap : atomicPropositions) {
                if (this->model.getStateLabeling().getStateHasLabel(ap, representativeState)) {
                    newLabeling.addLabelToState(ap, blockIndex);
                }
            }
        } else {
            // Compute the outgoing transitions of the block.
            std::map<storm::storage::sparse::state_type, ValueType> blockProbability;
            if (this->options.getUsesEpsilon()) {
                auto blockIterator = block.begin();

                // Start block distribution over partition with the first state of the block
                auto currentState = *blockIterator;
                auto blockDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
                for (auto const& e : this->model.getTransitionMatrix().getRow(currentState)) {
                    blockDistribution.addProbability(blocksMapping.at(this->partition.getBlockOfElement(e.getColumn()).front()), e.getValue());
                }
                ++blockIterator;

                // Now we enhance the block distribution with the intervals of every other state in the block
                while (blockIterator != block.end()) {
                    currentState = *blockIterator;
                    auto distributionOfState = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
                    for (auto const& e : this->model.getTransitionMatrix().getRow(currentState)) {
                        distributionOfState.addProbability(blocksMapping.at(this->partition.getBlockOfElement(e.getColumn()).front()), e.getValue());
                    }
                    blockDistribution = computeEnhancedDistribution(blockDistribution, distributionOfState);
                    ++blockIterator;
                }

                // Now add them to the actual matrix
                auto blockDistributionIterator = blockDistribution.begin();
                while (blockDistributionIterator != blockDistribution.end()) {
                    builder.addNextValue(blockIndex, blockDistributionIterator->first, blockDistributionIterator->second);
                    blockDistributionIterator++;
                }
            } else {
                // TODO: We could also use the feasible interval here for quotient construction.
                for (auto const& entry : this->model.getTransitionMatrix().getRow(representativeState)) {
                    auto targetBlock = blocksMapping.at(this->partition.getBlockOfElement(entry.getColumn()).front());

                    // TODO: Extend for weak bisimulation

                    auto probIterator = blockProbability.find(targetBlock);
                    if (probIterator != blockProbability.end()) {
                        probIterator->second += entry.getValue();
                        probIterator->second = clampIntervalToProbabilisticInterval(probIterator->second);
                    } else {
                        blockProbability[targetBlock] = entry.getValue();
                    }
                }

                // Now add them to the actual matrix
                for (auto const& probabilityEntry : blockProbability) {
                    // TODO: Handle case for weak bisimulation
                    builder.addNextValue(blockIndex, probabilityEntry.first, probabilityEntry.second);
                }
            }

            // Otherwise add all atomic propositions to the equivalence class that the representative state
            // satisfies.
            for (auto const& ap : atomicPropositions) {
                if (this->model.getStateLabeling().getStateHasLabel(ap, representativeState)) {
                    newLabeling.addLabelToState(ap, blockIndex);
                }
            }
        }

        // If the model has state rewards, we simply copy the state reward of the representative state, because
        // all states in a block are guaranteed to have the same state reward.
        if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
            auto const& rewardModel = this->model.getUniqueRewardModel();
            if (rewardModel.hasStateRewards()) {
                stateRewards.value()[blockIndex] = rewardModel.getStateRewardVector()[representativeState];
            }
            if (rewardModel.hasStateActionRewards()) {
                stateRewards.value()[blockIndex] += rewardModel.getStateActionRewardVector()[representativeState];
            }
        }
    });

    // Now check which of the blocks of the partition contain at least one initial state.
    for (auto initialState : this->model.getInitialStates()) {
        auto initialBlock = this->partition.getBlockOfElement(initialState);
        newLabeling.addLabelToState("init", blocksMapping.at(initialBlock.front()));
    }

    // Construct the reward model mapping.
    std::unordered_map<std::string, typename ModelType::RewardModelType> rewardModels;
    if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
        STORM_LOG_THROW(this->model.hasUniqueRewardModel(), storm::exceptions::IllegalFunctionCallException, "Cannot preserve more than one reward model.");
        auto nameRewardModelPair = this->model.getRewardModels().begin();
        rewardModels.insert(std::make_pair(nameRewardModelPair->first, typename ModelType::RewardModelType(stateRewards)));
    }

    // Finally construct the quotient model.
    this->quotient = std::make_shared<ModelType>(builder.build(), std::move(newLabeling), std::move(rewardModels));
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::postProcessInitialPartition() {
    for (int i = 0; i < this->model.getNumberOfStates(); i++) {
        auto distribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();

        for (auto const& e : this->model.getTransitionMatrix().getRow(i)) {
            distribution.addProbability(this->partition.getBlockOfElement(e.getColumn()).front(), e.getValue());
        }

        originalStateDistributions[i] = getClampedDistribution(distribution);
    }
}

template<typename ModelType>
storm::storage::Distribution<typename ModelType::ValueType, storm::storage::sparse::state_type>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::computeEnhancedDistribution(
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& firstDistribution,
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& secondDistribution) {
    auto enhancedDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
    std::set<storm::storage::sparse::state_type> keys;

    for (auto const& entry : firstDistribution) {
        keys.insert(entry.first);
    }

    for (auto const& entry : secondDistribution) {
        keys.insert(entry.first);
    }

    for (auto key : keys) {
        auto firstInterval = firstDistribution.getProbability(key);
        auto secondInterval = secondDistribution.getProbability(key);

        auto lower = firstInterval.lower() < secondInterval.lower() ? firstInterval.lower() : secondInterval.lower();
        auto upper = firstInterval.upper() > secondInterval.upper() ? firstInterval.upper() : secondInterval.upper();

        enhancedDistribution.addProbability(key, ValueType(lower, carl::BoundType::WEAK, upper, carl::BoundType::WEAK));
    }

    return getClampedDistribution(enhancedDistribution);
}

template<typename ModelType>
storm::IntervalBaseType<typename ModelType::ValueType> DeterministicIntervalModelBisimulationDecomposition<ModelType>::computeDeltaForState(
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& stateDistribution,
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& enhancedDistribution) {
    storm::IntervalBaseType<ValueType> delta = 0.0;

    auto enhancedDistributionIterator = enhancedDistribution.begin();
    while (enhancedDistributionIterator != enhancedDistribution.end()) {
        ValueType intervalFromState = stateDistribution.getProbability(enhancedDistributionIterator->first);
        ValueType intervalFromEnhancedState = enhancedDistributionIterator->second;

        IntervalBaseType<ValueType> deltaOfLowerBound = intervalFromState.lower() - intervalFromEnhancedState.lower();
        IntervalBaseType<ValueType> deltaOfUpperBound = intervalFromState.upper() - intervalFromEnhancedState.upper();
        delta += (storm::utility::abs(deltaOfLowerBound)) + (storm::utility::abs(deltaOfUpperBound));

        enhancedDistributionIterator++;
    }

    return delta;
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::debugFindMergeableFinalBlocks(double epsilon) {
    std::vector<storm::storage::Partition::Block> blocks;
    this->partition.forEachBlock([&](auto const& block) { blocks.push_back(block); });

    for (std::size_t i = 0; i < blocks.size(); ++i) {
        for (std::size_t j = i + 1; j < blocks.size(); ++j) {
            if (blocks[i].empty() || blocks[j].empty()) {
                continue;
            }

            auto si = blocks[i].front();
            auto sj = blocks[j].front();

            bool sameLabels = true;
            for (auto const& ap : this->options.respectedAtomicPropositions.value()) {
                if (this->model.getStateLabeling().getStateHasLabel(ap, si) != this->model.getStateLabeling().getStateHasLabel(ap, sj)) {
                    sameLabels = false;
                    break;
                }
            }

            if (!sameLabels) {
                continue;
            }

            if (canMergeBlocksForDebug(blocks[i], blocks[j], epsilon)) {
                std::cout << "DEBUG mergeable final blocks: " << blocks[i].front() << " and " << blocks[j].front() << " with sizes " << blocks[i].size()
                          << " and " << blocks[j].size() << std::endl;
            }
        }
    }
}

template<typename ModelType>
bool DeterministicIntervalModelBisimulationDecomposition<ModelType>::canMergeBlocksForDebug(std::span<uint64_t const> blockA, std::span<uint64_t const> blockB,
                                                                                            double epsilon) {
    std::vector<storm::storage::sparse::state_type> states;
    states.insert(states.end(), blockA.begin(), blockA.end());
    states.insert(states.end(), blockB.begin(), blockB.end());

    if (states.empty()) {
        return true;
    }

    auto candidate = originalStateDistributions[states.front()];

    for (std::size_t i = 1; i < states.size(); ++i) {
        candidate = computeEnhancedDistribution(candidate, originalStateDistributions[states[i]]);
    }

    bool ok = true;
    for (auto s : states) {
        auto delta = computeDeltaForState(originalStateDistributions[s], candidate);

        // std::cout << "merge-debug state " << s << ", delta = " << delta << ", 2delta = " << (2 * delta) << ", epsilon = " << epsilon << std::endl;

        if (2 * delta > epsilon) {
            ok = false;
        }
    }

    return ok;
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::validateEpsilonStability(double epsilon) {
    bool valid = true;

    this->partition.forEachBlock([&](auto const& block) {
        if (block.empty() || this->absorbingBlocks.contains(block.front())) {
            return;
        }

        auto candidate = originalStateDistributions[block.front()];

        for (auto it = std::next(block.begin()); it != block.end(); ++it) {
            candidate = computeEnhancedDistribution(candidate, originalStateDistributions[*it]);
        }

        for (auto state : block) {
            auto delta = computeDeltaForState(originalStateDistributions[state], candidate);

            if (2 * delta > epsilon) {
                valid = false;
                std::cout << "EPSILON VALIDATION FAILED: block front=" << block.front() << ", state=" << state << ", delta=" << delta
                          << ", 2delta=" << (2 * delta) << ", epsilon=" << epsilon << std::endl;
            }
        }
    });

    if (valid) {
        std::cout << "Epsilon validation passed for final partition." << std::endl;
    }
}

template<typename ModelType>
typename ModelType::ValueType DeterministicIntervalModelBisimulationDecomposition<ModelType>::clampIntervalToProbabilisticInterval(ValueType interval) const {
    auto lowerBound = interval.lower();
    auto upperBound = interval.upper();

    auto const zero = storm::utility::zero<IntervalBaseType<ValueType>>();
    auto const one = storm::utility::one<IntervalBaseType<ValueType>>();

    if (!this->model.isExact()) {
        auto tolerance = this->options.getTolerance();

        if (tolerance.isPointInterval()) {
            if (lowerBound > one && lowerBound - one <= tolerance.lower()) {
                lowerBound = one;
            }

            if (upperBound > one && upperBound - one <= tolerance.lower()) {
                upperBound = one;
            }

            if (lowerBound > upperBound) {
                STORM_LOG_ERROR("Applying tolerance led to an invalid interval!");
            }

            if (lowerBound > 1 || upperBound > 1) {
                STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Computation led to an invalid interval!");
            }
        } else {
            STORM_LOG_ERROR("Unable to apply interval tolerance to non-exact model, as tolerance is not a point-interval!");
        }
    }

    auto zeroOneInterval = ValueType(zero, carl::BoundType::WEAK, one, carl::BoundType::WEAK);

    return ValueType(lowerBound, carl::BoundType::WEAK, upperBound, carl::BoundType::WEAK).intersect(zeroOneInterval);
}

template<typename ModelType>
bool DeterministicIntervalModelBisimulationDecomposition<ModelType>::checkCurrentPartitionByExactFeasibleIntervals() const {
    auto computeFeasibleInterval = [&](storm::storage::sparse::state_type state, auto const& splitterBlock) {
        ValueType toBlock = storm::utility::zero<ValueType>();
        ValueType toOther = storm::utility::zero<ValueType>();

        for (auto const& entry : this->model.getTransitionMatrix().getRow(state)) {
            if (this->partition.contains(splitterBlock, entry.getColumn())) {
                toBlock += entry.getValue();
            } else {
                toOther += entry.getValue();
            }
        }

        return computeFeasibleIntervalBasedOnAggregatedIntervals(clampIntervalToProbabilisticInterval(toBlock), clampIntervalToProbabilisticInterval(toOther));
    };

    bool valid = true;

    this->partition.forEachBlock([&](auto const& block) {
        if (!valid || block.size() <= 1) {
            return;
        }

        auto representative = block.front();

        for (auto state : block) {
            // Check labels.
            for (auto const& ap : this->options.respectedAtomicPropositions.value()) {
                bool repHasLabel = this->model.getStateLabeling().getStateHasLabel(ap, representative);
                bool stateHasLabel = this->model.getStateLabeling().getStateHasLabel(ap, state);

                if (repHasLabel != stateHasLabel) {
                    valid = false;
                    std::cout << "Exact feasible-interval partition check failed: labels differ in block " << block.front()
                              << ", representative=" << representative << ", state=" << state << ", AP=" << ap << std::endl;
                    return;
                }
            }

            // Check feasible intervals against every splitter block.
            this->partition.forEachBlock([&](auto const& splitterBlock) {
                if (!valid) {
                    return;
                }

                auto repInterval = computeFeasibleInterval(representative, splitterBlock);
                auto stateInterval = computeFeasibleInterval(state, splitterBlock);

                if (repInterval != stateInterval) {
                    valid = false;

                    std::cout << "Exact feasible-interval partition check failed." << std::endl;
                    std::cout << "  block representative: " << representative << std::endl;
                    std::cout << "  state: " << state << std::endl;
                    std::cout << "  splitter block representative: " << splitterBlock.front() << std::endl;
                    std::cout << "  representative interval: " << repInterval << std::endl;
                    std::cout << "  representative distribution: " << this->originalStateDistributions[representative] << std::endl;
                    std::cout << "  state distribution: " << this->originalStateDistributions[state] << std::endl;
                    std::cout << "  state interval: " << stateInterval << std::endl;
                }
            });

            if (!valid) {
                return;
            }
        }
    });

    if (valid) {
        std::cout << "Exact feasible-interval partition check passed." << std::endl;
    }

    return valid;
}

template<typename ModelType>
std::pair<storm::storage::BitVector, storm::storage::BitVector> DeterministicIntervalModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    // TODO: return actual state sets for probability 0 and 1
    return {storm::storage::BitVector(), storm::storage::BitVector()};
}

template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;
template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::Interval>>;
template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>;
template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalInterval>>;

}  // namespace storage
}  // namespace storm

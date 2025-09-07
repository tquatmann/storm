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
      touchedProbabilitiesToSplitter(model.getNumberOfStates(), false) {}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::buildQuotient() {
    std::cout << "Size of partition: " << this->partition.getNumberOfBlocks() << std::endl;

    // In order to create the quotient model, we need to construct
    // (a) the new transition matrix,
    // (b) the new labeling,
    // (c) the new reward structures.

    // Prepare a matrix builder for (a).
    storm::storage::SparseMatrixBuilder<ValueType> builder(this->partition.getNumberOfBlocks(), this->partition.getNumberOfBlocks());

    // Prepare the new state labeling for (b).
    storm::models::sparse::StateLabeling newLabeling(this->partition.getNumberOfBlocks());
    std::set<std::string> atomicPropositionsSet = this->options.respectedAtomicPropositions.get();
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
            for (auto const& entry : this->model.getTransitionMatrix().getRow(representativeState)) {
                auto targetBlock = blocksMapping.at(this->partition.getBlockOfElement(entry.getColumn()).front());

                // If we are computing a weak bisimulation quotient, there is no need to add self-loops.
                // TODO: Extend for weak bisimulation

                auto probIterator = blockProbability.find(targetBlock);
                if (probIterator != blockProbability.end()) {
                    probIterator->second += entry.getValue();

                    // Normalize interval
                    auto normalizedInterval = storm::Interval(std::min(1.0, probIterator->second.lower()), std::min(1.0, probIterator->second.upper()));
                    probIterator->second = normalizedInterval;
                } else {
                    blockProbability[targetBlock] = entry.getValue();
                }
            }

            // Now add them to the actual matrix.
            for (auto const& probabilityEntry : blockProbability) {
                // TODO: Handle case for weak bisimulation
                builder.addNextValue(blockIndex, probabilityEntry.first, probabilityEntry.second);
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
    this->quotient = std::shared_ptr<ModelType>(new ModelType(builder.build(), std::move(newLabeling), std::move(rewardModels)));
}

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
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::postProcessInitialPartition() {
    // TODO: implement
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(
    storm::storage::Partition::Block splitterBlock,
    std::deque<typename bisimulation::Partition::Block>& splitterQueue,
    bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::bisimulation::Partition::BlockSet blocksToSplit;
    // std::cout << "Performing interval bisimulation!" << std::endl;

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
                probabilitiesToOtherBlocks[predecessorState] = storm::utility::zero<ValueType>(); // init I(s, \Pi \setminus C)
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
                probabilitiesToOtherBlocks[predecessorState] += entry.getValue(); // interval addition
            }
        }
    }

    for (auto predecessorBlockToSplit : blocksToSplit) {
        // First split the block by whether it is a predecessor of the splitter block or not
        auto [noPredecessors, predecessors] =
            this->partition.splitBlockByPredicate(predecessorBlockToSplit, [this]
                                                  (auto const& state) { return touchedProbabilitiesToSplitter.get(state); });

        if (noPredecessors.size() > 0) {
            if (!enqueuedSplitterBlocks.contains(noPredecessors)) {
                splitterQueue.push_back(noPredecessors);
                enqueuedSplitterBlocks.insert(noPredecessors);
            }
        }

        if (predecessors.size() > 0) {
            bool wasSplit = this->partition.splitBlockByOrder(predecessors, [this]
                                                              (auto const& a, auto const& b) {
                                                                  if (!touchedProbabilitiesToSplitter.get(a) || !touchedProbabilitiesToSplitter.get(b)) {
                                                                      std::cout << "Comparing states that were not touched: " << a << ", " << b << std::endl;
                                                                      return false;
                                                                  }

                                                                  auto projectedIntervalOfStateA = computeIntervalProjection(probabilitiesToCurrentSplitter[a], probabilitiesToOtherBlocks[a]);
                                                                  auto projectedIntervalOfStateB = computeIntervalProjection(probabilitiesToCurrentSplitter[b], probabilitiesToOtherBlocks[b]);

                                                                  // TODO: [0.1, 0.2] < [0.1, 0.3] returns false!
                                                                  // auto result = projectedIntervalOfStateA < projectedIntervalOfStateB;

                                                                  return projectedIntervalOfStateA.lower() < projectedIntervalOfStateB.lower() ||
                                                                         (projectedIntervalOfStateA.lower() == projectedIntervalOfStateB.lower() &&
                                                                          projectedIntervalOfStateA.upper() < projectedIntervalOfStateB.upper());
                                                              });

            // Add all blocks that were split to splitter queue
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

    // Reset the touched entries of the probabilitiesToCurrentSplitter vector
    touchedProbabilitiesToSplitter.clear();
}


template<typename ModelType>
ModelType::ValueType DeterministicIntervalModelBisimulationDecomposition<ModelType>::computeIntervalProjection(ValueType intervalToSplitter, ValueType intervalToOtherBlocks) {
    // Normalize intervals
    carl::Interval normalizedIntervalToSplitter(std::min(1.0, intervalToSplitter.lower()), std::min(1.0, intervalToSplitter.upper()));
    carl::Interval normalizedIntervalToOtherBlocks(std::min(1.0, intervalToOtherBlocks.lower()), std::min(1.0, intervalToOtherBlocks.upper()));

    // Compute interval projection
    return carl::Interval(std::max(normalizedIntervalToSplitter.lower(), 1.0 - normalizedIntervalToOtherBlocks.upper()),
                          std::min(normalizedIntervalToSplitter.upper(), 1.0 - normalizedIntervalToOtherBlocks.lower()));
}

template<typename ModelType>
bool DeterministicIntervalModelBisimulationDecomposition<ModelType>::possiblyNeedsRefinement(std::span<uint64_t const> block) const {
    return block.size() > 1 && !this->absorbingBlocks.contains(block.front());
}

template<typename ModelType>
std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>>
    DeterministicIntervalModelBisimulationDecomposition<ModelType>::create2DPolytope(storm::RationalNumber c1LowerBound, storm::RationalNumber c1UpperBound,
                                                                                            storm::RationalNumber c2LowerBound, storm::RationalNumber c2UpperBound) {
    using Point = typename storm::storage::geometry::Polytope<storm::RationalNumber>::Point;

    // Halfspace in R² is given by: a1 * p1 + a2 * p2 <= b, where (a1, a2) is the normal vector and b is the offset
    std::vector<storm::storage::geometry::Halfspace<storm::RationalNumber>> halfspaces;

    // Interval bounds
    halfspaces.emplace_back(Point{-1.0, 0.0}, // a1, a2
                            -c1LowerBound);  // => (-p1 <= -l1) <=> (p1 >= l1)
    halfspaces.emplace_back(Point{ 1.0, 0.0}, // a1, a2
                            c1UpperBound);  // => (p1 <= u1)
    halfspaces.emplace_back(Point{ 0.0,-1.0}, // a1, a2
                            -c2LowerBound);  // => (-p2 <= -l2) <=> (p2 >= l2)
    halfspaces.emplace_back(Point{ 0.0, 1.0}, // a1, a2
                            c2UpperBound);  // => (p2 <= u2)

    // Normalization constraint: p1 + p2 = 1
    halfspaces.emplace_back(Point{ 1.0, 1.0}, // a1, a2
                            1.0);  // p1 + p2 <= 1
    halfspaces.emplace_back(Point{-1.0,-1.0}, // a1, a2
                            -1.0);  // -(p1 + p2) <= -1  => p1 + p2 >= 1

    // Polytope is "intersection" of all halfspaces
    return storm::storage::geometry::Polytope<storm::RationalNumber>::create(halfspaces);
}

template<typename ModelType>
std::pair<storm::storage::BitVector, storm::storage::BitVector>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    // TODO: return actual state sets for probability 0 and 1
    return {storm::storage::BitVector(), storm::storage::BitVector()};
}

template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<carl::Interval<double>>>;
template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Ctmc<carl::Interval<double>>>;

}  // namespace storage
}  // namespace storm

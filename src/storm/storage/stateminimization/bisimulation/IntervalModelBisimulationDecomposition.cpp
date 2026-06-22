#include "IntervalModelBisimulationDecomposition.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
IntervalModelBisimulationDecomposition<ModelType>::IntervalModelBisimulationDecomposition(
    const ModelType& model, const typename BisimulationDecomposition<ModelType>::BisimulationOptions& options)
    : storm::storage::BisimulationDecomposition<ModelType>(model, options),
      probabilitiesToCurrentSplitter(model.getTransitionMatrix().getRowCount(), storm::utility::zero<ValueType>()),
      probabilitiesToOtherBlocks(model.getTransitionMatrix().getRowCount(), storm::utility::zero<ValueType>()),
      touchedProbabilitiesToSplitter(model.getTransitionMatrix().getRowCount(), false),
      choiceToStateMapping(model.getTransitionMatrix().getRowCount()) {
    // Initialize choice to state mapping to compute interval probabilities to the splitter blocks later.
    auto const& choiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    for (auto state = 0; state < this->model.getNumberOfStates(); ++state) {
        for (auto choice = choiceIndices[state]; choice < choiceIndices[state + 1]; ++choice) {
            choiceToStateMapping[choice] = state;
        }
    }

    backwardTransitionsWithChoices = this->model.getTransitionMatrix().transpose(false);
}

template<typename ModelType>
void IntervalModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(storm::storage::stateminimization::Partition::Block splitterBlock,
                                                                                       std::deque<typename stateminimization::Partition::Block>& splitterQueue,
                                                                                       stateminimization::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::stateminimization::Partition::BlockSet blocksToSplit;
    auto choiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();

    for (auto currentState : splitterBlock) {
        // Compute the transition intervals for each choice/action.
        // Note that for IDTMCs the choice distribution entry equals the state distribution entry.
        for (const auto& predecessorChoiceEntry : backwardTransitionsWithChoices.getRow(currentState)) {
            auto predecessorChoice = predecessorChoiceEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(choiceToStateMapping[predecessorChoice]);

            auto intervalTransitionProbability = predecessorChoiceEntry.getValue();

            if (!possiblyNeedsRefinement(predecessorBlock)) {
                continue;
            }

            // Compute interval going to the splitter, i.e., I(s, C), and initialize interval leading to all other blocks, i.e., I(s, \Pi \setminus C).
            if (touchedProbabilitiesToSplitter.get(predecessorChoice)) {
                probabilitiesToCurrentSplitter[predecessorChoice] += intervalTransitionProbability;
            } else {
                probabilitiesToCurrentSplitter[predecessorChoice] = intervalTransitionProbability;
                probabilitiesToOtherBlocks[predecessorChoice] = storm::utility::zero<ValueType>();
                touchedProbabilitiesToSplitter.set(predecessorChoice, true);
            }

            // Remember which blocks contain predecessors to split w.r.t. the splitter afterwards
            blocksToSplit.insert(predecessorBlock);
        }
    }

    // Compute interval to all blocks except the splitter, i.e., I(s, \alpha, \Pi \setminus C) for IMDPs or I(s, \Pi \setminus C) for IDTMCs.
    for (auto predecessorChoice : touchedProbabilitiesToSplitter) {
        for (const auto& entry : this->model.getTransitionMatrix().getRow(predecessorChoice)) {
            auto targetState = entry.getColumn();

            if (!this->partition.contains(splitterBlock, targetState)) {
                probabilitiesToOtherBlocks[predecessorChoice] += entry.getValue();
            }
        }
    }

    for (auto predecessorBlockToSplit : blocksToSplit) {
        // First split the block by whether it is a predecessor of the splitter block or not.
        auto [noPredecessors, predecessors] = this->partition.splitBlockByPredicate(predecessorBlockToSplit, [this, &choiceIndices](auto const& state) {
            for (auto choice = choiceIndices[state]; choice < choiceIndices[state + 1]; choice++) {
                // Check if there is any choice of a state that leads to a state in the splitter.
                if (touchedProbabilitiesToSplitter.get(choice)) {
                    return true;
                }
            }
            return false;
        });

        if (noPredecessors.size() > 0) {
            if (!enqueuedSplitterBlocks.contains(noPredecessors)) {
                splitterQueue.push_back(noPredecessors);
                enqueuedSplitterBlocks.insert(noPredecessors);
            }
        }

        if (predecessors.size() > 0) {
            bool wasSplit = this->partition.splitBlockByOrder(predecessors, [this, &choiceIndices](auto const& a, auto const& b) {
                auto firstChoiceOfStateA = choiceIndices[a];
                auto lastChoiceOfStateA = choiceIndices[a + 1];
                auto firstChoiceOfStateB = choiceIndices[b];
                auto lastChoiceOfStateB = choiceIndices[b + 1];

                STORM_LOG_ASSERT(lastChoiceOfStateA - firstChoiceOfStateA == lastChoiceOfStateB - firstChoiceOfStateB,
                                 "States with different numbers of choices should have been split before interval refinement.");

                for (auto choiceOffset = 0; choiceOffset < lastChoiceOfStateA - firstChoiceOfStateA; ++choiceOffset) {
                    auto currentChoiceOfStateA = firstChoiceOfStateA + choiceOffset;
                    auto currentChoiceOfStateB = firstChoiceOfStateB + choiceOffset;

                    auto feasibleIntervalOfChoiceA = computeFeasibleIntervalBasedOnAggregatedIntervals(probabilitiesToCurrentSplitter[currentChoiceOfStateA],
                                                                                                       probabilitiesToOtherBlocks[currentChoiceOfStateA]);
                    auto feasibleIntervalOfChoiceB = computeFeasibleIntervalBasedOnAggregatedIntervals(probabilitiesToCurrentSplitter[currentChoiceOfStateB],
                                                                                                       probabilitiesToOtherBlocks[currentChoiceOfStateB]);

                    if (!this->comparator.isEqual(feasibleIntervalOfChoiceA.lower(), feasibleIntervalOfChoiceB.lower())) {
                        return this->comparator.isLess(feasibleIntervalOfChoiceA.lower(), feasibleIntervalOfChoiceB.lower());
                    }
                    if (!this->comparator.isEqual(feasibleIntervalOfChoiceA.upper(), feasibleIntervalOfChoiceB.upper())) {
                        return this->comparator.isLess(feasibleIntervalOfChoiceA.upper(), feasibleIntervalOfChoiceB.upper());
                    }
                }

                return false;
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
ModelType::ValueType IntervalModelBisimulationDecomposition<ModelType>::computeFeasibleIntervalBasedOnAggregatedIntervals(
    ValueType intervalToSplitter, ValueType intervalToOtherBlocks) const {
    // Clamp intervals.
    ValueType clampedIntervalToSplitter = storm::utility::interval::makeIntervalProbabilistic(intervalToSplitter, this->comparator);
    ValueType clampedIntervalToOtherBlocks = storm::utility::interval::makeIntervalProbabilistic(intervalToOtherBlocks, this->comparator);

    // Compute feasible interval.
    storm::IntervalBaseType<ValueType> lowerBoundOfFeasibleInterval;
    if (this->comparator.isLess(storm::utility::one<IntervalBaseType<ValueType>>() - clampedIntervalToOtherBlocks.upper(), clampedIntervalToSplitter.lower())) {
        lowerBoundOfFeasibleInterval = clampedIntervalToSplitter.lower();
    } else {
        lowerBoundOfFeasibleInterval = storm::utility::one<IntervalBaseType<ValueType>>() - clampedIntervalToOtherBlocks.upper();
    }

    storm::IntervalBaseType<ValueType> upperBoundOfFeasibleInterval;
    if (this->comparator.isLess(clampedIntervalToSplitter.upper(), storm::utility::one<IntervalBaseType<ValueType>>() - clampedIntervalToOtherBlocks.lower())) {
        upperBoundOfFeasibleInterval = clampedIntervalToSplitter.upper();
    } else {
        upperBoundOfFeasibleInterval = storm::utility::one<IntervalBaseType<ValueType>>() - clampedIntervalToOtherBlocks.lower();
    }

    // For non-exact computations, it might be that the lower bound is slightly larger than the upper bound due to imprecision.
    // Thus, we make sure here to make them equal to avoid returning an empty interval.
    if (this->comparator.isEqual(lowerBoundOfFeasibleInterval, upperBoundOfFeasibleInterval)) {
        upperBoundOfFeasibleInterval = lowerBoundOfFeasibleInterval;
    }

    return ValueType(lowerBoundOfFeasibleInterval, carl::BoundType::WEAK, upperBoundOfFeasibleInterval, carl::BoundType::WEAK);
}

template<typename ModelType>
bool IntervalModelBisimulationDecomposition<ModelType>::possiblyNeedsRefinement(std::span<uint64_t const> block) const {
    return block.size() > 1 && !this->absorbingBlocks.contains(block.front());
}

template<typename ModelType>
void IntervalModelBisimulationDecomposition<ModelType>::buildQuotientFromPartition() {
    // In order to create the quotient model, we need to construct
    // (a) the new transition matrix,
    // (b) the new labeling,
    // (c) the new reward structures.

    // Prepare a matrix builder for (a).
    storm::storage::SparseMatrixBuilder<ValueType> builder(0, this->partition.getNumberOfBlocks(), 0, false, true, this->partition.getNumberOfBlocks());

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

    // If the model had action rewards, we need to build the state rewards for the quotient as well.
    std::optional<std::vector<ValueType>> stateActionRewards;
    if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
        if (this->model.getUniqueRewardModel().hasStateActionRewards()) {
            stateActionRewards = std::vector<ValueType>();
        }
    }

    // Now build (a) and (b) by traversing all blocks.
    std::map<uint64_t, uint64_t> blocksMapping;
    auto blockIndex = 0;
    this->partition.forEachBlock([&](auto const& block) {
        blocksMapping.emplace(block.front(), blockIndex);
        blockIndex++;
    });

    uint_fast64_t currentRow = 0;
    this->partition.forEachBlock([&](auto const& block) {
        builder.newRowGroup(currentRow);

        // Pick representative state.
        // Here, it does not matter which one to take, as we are constructing the enhanced intervals of the quotient state based on all states of the block.
        auto representativeState = block.front();
        blockIndex = blocksMapping.at(representativeState);

        // If the block is absorbing, we simply add a self-loop.
        if (this->absorbingBlocks.contains(representativeState)) {
            builder.addNextValue(blockIndex, blockIndex, storm::utility::one<ValueType>());
            ++currentRow;

            // If the block has a special representative state, we retrieve it now.
            representativeState = this->absorbingBlocks.at(representativeState);

            // Give the choice a reward of zero as we artificially introduced that the block is absorbing.
            if (this->options.getKeepRewards() && this->model.hasRewardModel() && this->model.getUniqueRewardModel().hasStateActionRewards()) {
                stateActionRewards.value().push_back(storm::utility::zero<ValueType>());
            }

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
            auto choiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();
            auto numberOfChoicesInBlock = choiceIndices[representativeState + 1] - choiceIndices[representativeState];
            for (uint_fast64_t choiceOffset = 0; choiceOffset < numberOfChoicesInBlock; choiceOffset++) {
                // Start block distribution for choice over final partition with the first state of the block.
                auto blockStateIterator = block.begin();
                auto currentBlockState = *blockStateIterator;
                auto blockChoiceDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();

                // TODO: Do we want to have the feasible intervals for the quotient?
                // Compute distribution of 'currentBlockState' over the quotient partition with current choice.
                for (auto const& e : this->model.getTransitionMatrix().getRow(choiceIndices[currentBlockState] + choiceOffset)) {
                    blockChoiceDistribution.addProbability(blocksMapping.at(this->partition.getBlockOfElement(e.getColumn()).front()), e.getValue());
                }

                blockChoiceDistribution = storm::utility::interval::makeDistributionProbabilistic(blockChoiceDistribution, this->comparator);

                // Now add them to the actual matrix.
                auto blockDistributionIterator = blockChoiceDistribution.begin();
                while (blockDistributionIterator != blockChoiceDistribution.end()) {
                    builder.addNextValue(currentRow, blockDistributionIterator->first, blockDistributionIterator->second);
                    blockDistributionIterator++;
                }

                // TODO: Implement this by using storm::storage::DistributionWithReward instead of storm::storage::Distribution.
                // if (this->options.getKeepRewards() && this->model.hasRewardModel() && this->model.getUniqueRewardModel().hasStateActionRewards()) {
                //     stateActionRewards.value().push_back(quotientDistributions[choice].getReward());
                // }
                ++currentRow;
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

    // In case of deterministic model, make row grouping trivial.
    auto transitionMatrix = builder.build();
    if (!this->model.isNondeterministicModel()) {
        transitionMatrix.makeRowGroupingTrivial();
    }

    // Finally construct the quotient model.
    this->quotient = std::make_shared<ModelType>(std::move(transitionMatrix), std::move(newLabeling), std::move(rewardModels));
}

// TODO: Check this in an assert or during tests?
// TODO: Need to adjust this first though.
template<typename ModelType>
bool IntervalModelBisimulationDecomposition<ModelType>::checkCurrentPartitionByExactFeasibleIntervals() const {
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

        return computeFeasibleIntervalBasedOnAggregatedIntervals(storm::utility::interval::makeIntervalProbabilistic(toBlock),
                                                                 storm::utility::interval::makeIntervalProbabilistic(toOther));
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
                    // std::cout << "  representative distribution: " << this->originalChoiceDistributions[representative] << std::endl;
                    // std::cout << "  state distribution: " << this->originalChoiceDistributions[state] << std::endl;
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
std::pair<storm::storage::BitVector, storm::storage::BitVector> IntervalModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    // TODO: return actual state sets for probability 0 and 1
    return {storm::storage::BitVector(), storm::storage::BitVector()};
}

template class IntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;
template class IntervalModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::Interval>>;
template class IntervalModelBisimulationDecomposition<storm::models::sparse::Mdp<storm::Interval>>;
template class IntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>;
template class IntervalModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalInterval>>;
template class IntervalModelBisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalInterval>>;

}  // namespace storage
}  // namespace storm

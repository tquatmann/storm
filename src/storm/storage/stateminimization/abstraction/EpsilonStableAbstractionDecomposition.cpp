#include "storm/storage/stateminimization/abstraction/EpsilonStableAbstractionDecomposition.h"

namespace storm {
namespace storage {
namespace abstraction {

using namespace abstraction;

template<typename ModelType>
EpsilonStableAbstractionDecomposition<ModelType>::EpsilonStableAbstractionOptions::EpsilonStableAbstractionOptions(double epsilon)
    : BaseDecomposition<ModelType>::BaseOptions(), epsilon(epsilon) {
    // Intentionally left empty.
}

template<typename ModelType>
EpsilonStableAbstractionDecomposition<ModelType>::EpsilonStableAbstractionOptions::EpsilonStableAbstractionOptions(ModelType const& model,
                                                                                                                   storm::logic::Formula const& formula,
                                                                                                                   double epsilon)
    : BaseDecomposition<ModelType>::BaseOptions(model, formula), epsilon(epsilon) {
    // Intentionally left empty.
}

template<typename ModelType>
EpsilonStableAbstractionDecomposition<ModelType>::EpsilonStableAbstractionOptions::EpsilonStableAbstractionOptions(
    ModelType const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas, double epsilon)
    : BaseDecomposition<ModelType>::BaseOptions(model, formulas), epsilon(epsilon) {
    // Intentionally left empty.
}

template<typename ModelType>
EpsilonStableAbstractionDecomposition<ModelType>::EpsilonStableAbstractionDecomposition(
    const ModelType& model, const EpsilonStableAbstractionDecomposition::EpsilonStableAbstractionOptions& options)
    : BaseDecomposition<ModelType>(model, model.getBackwardTransitions(), options.getTolerance()),
      options(options),
      originalChoiceDistributions(model.getTransitionMatrix().getRowCount()) {}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::computeInitialPartition() {
    this->splitInitialPartitionBasedOnLabelsAndRewards();
    if (this->model.isNondeterministicModel()) {
        this->splitInitialPartitionBasedOnActionSets();
    }
    this->initializeChoiceDistributions();
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::refinePartition() {
    this->performEpsilonStableAbstractionRefinement(this->options.getEpsilon());
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::performEpsilonStableAbstractionRefinement(double epsilon) {
    std::deque<typename stateminimization::Partition::Block> blocksQueue;
    storm::storage::stateminimization::Partition::BlockSet enqueuedBlocks;

    // Enqueue all non-absorbing blocks for refinement.
    this->partition.forEachBlock([&](auto const& block) {
        if (!this->absorbingBlocks.contains(block.front())) {
            blocksQueue.push_back(block);
            enqueuedBlocks.insert(block);
        }
    });

    // Refinement loop
    uint_fast64_t iterations = 0;
    while (!blocksQueue.empty()) {
        ++iterations;

        auto block = blocksQueue.back();
        blocksQueue.pop_back();
        enqueuedBlocks.erase(block);

        refineBlockBasedOnEpsilonSignature(block, blocksQueue, enqueuedBlocks, epsilon);
    }
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::refineBlockBasedOnEpsilonSignature(std::span<uint64_t const> block,
                                                                                          std::deque<typename stateminimization::Partition::Block>& blocksQueue,
                                                                                          stateminimization::Partition::BlockSet& enqueuedBlocks,
                                                                                          double epsilon) {
    // Cluster states of this block in groups.
    std::vector<EpsilonStableAbstractionDecomposition::StateGroup> groups;

    // TODO: Can we assume that choices/actions, that are labeled the same, are actually stored in the same order in the rowGrouping for each state?
    for (auto s : block) {
        bool placed = false;
        auto distributionsOfS = getChoiceDistributionsOfState(s);
        for (auto& g : groups) {
            bool fits = true;

            STORM_LOG_ASSERT(distributionsOfS.size() == g.getNumberOfChoices(), "States with different number of choices/actions in one block.");

            std::vector<storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>> candidateDistributions;
            candidateDistributions.reserve(g.getNumberOfChoices());

            // Iterate over all possible choices/actions and check the \varepsilon-stable criterion for each.
            for (std::size_t choiceOffset = 0; choiceOffset < distributionsOfS.size(); ++choiceOffset) {
                // Compute resulting enhanced distribution of current choice/action on adding state 's' to the group 'g'.
                auto candidateDistribution = this->computeEnhancedDistribution(g.getDistributionAtOffset(choiceOffset), distributionsOfS[choiceOffset]);

                // Check if enhancement satisfies \varepsilon-criterion via our upper bound w.r.t. s itself
                if (this->comparator.isLess(epsilon, 2 * computeDeltaForState(distributionsOfS[choiceOffset], candidateDistribution))) {
                    // If we cannot satisfy our criterion for one choice/action, then the state does not fit in the current group 'g'.
                    fits = false;
                    break;
                }

                for (auto t : g.getStates()) {
                    auto distributionsOfT = getChoiceDistributionsOfState(t);

                    STORM_LOG_ASSERT(distributionsOfT.size() == distributionsOfS.size(),
                                     "States in an epsilon-stable group have incompatible numbers of choices.");

                    // Compute 'delta' between the state distribution of 't' and the candidate distribution of the current group distribution including the
                    // interval extensions based on candidate state 't'.
                    if (this->comparator.isLess(epsilon, 2 * computeDeltaForState(distributionsOfT[choiceOffset], candidateDistribution))) {
                        // If we cannot satisfy our criterion for one choice/action, then the state does not fit in the current group 'g'.
                        fits = false;
                        break;
                    }
                }

                // \varepsilon-criterion would be violated by adding 's' to group 'g'. Thus, we abort.
                if (!fits) {
                    break;
                }

                // Temporarily save candidate distribution for current choice/action.
                // We might need it later if the criterion is satisfied for each of the enabled choices/actions.
                candidateDistributions.push_back(std::move(candidateDistribution));
            }

            // Add state to current group 'g'
            if (fits) {
                g.addState(s);

                // Update enhanced distributions of group.
                for (std::size_t choiceOffset = 0; choiceOffset < candidateDistributions.size(); ++choiceOffset) {
                    g.setDistributionAtOffset(choiceOffset, std::move(candidateDistributions[choiceOffset]));
                }

                placed = true;
                break;
            }
        }

        // If the state was not placed yet, create its own group.
        if (!placed) {
            groups.push_back(EpsilonStableAbstractionDecomposition::StateGroup({s}, distributionsOfS));
        }
    }

    if (groups.size() <= 1)
        return;  // Nothing to split

    // Split block by grouping.
    // TODO: We might be able to implement this more efficiently, e.g., by using a stateToGroup mapping that gets filled while assigning the states to their
    // TODO: groups.
    bool wasSplit = this->partition.splitBlockByOrder(block, [&](auto a, auto b) {
        int groupOfA = -1, groupOfB = -1;
        for (int i = 0; i < groups.size(); ++i) {
            if (std::find(groups[i].getStates().begin(), groups[i].getStates().end(), a) != groups[i].getStates().end())
                groupOfA = i;
            if (std::find(groups[i].getStates().begin(), groups[i].getStates().end(), b) != groups[i].getStates().end())
                groupOfB = i;
        }
        return (groupOfA < groupOfB);
    });

    if (wasSplit) {
        // Place blocks of predecessors into queue.
        this->partition.forEachSubBlock(block, [this, &blocksQueue, &enqueuedBlocks](auto const& subBlock) {
            for (auto state : subBlock) {
                for (auto& transition : this->backwardTransitions.getRow(state)) {
                    auto predecessorState = transition.getColumn();
                    auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);

                    // Place predecessor block into queue only if it is not already there.
                    if (predecessorBlock.size() > 1 && !enqueuedBlocks.contains(predecessorBlock)) {
                        blocksQueue.push_back(predecessorBlock);
                        enqueuedBlocks.insert(predecessorBlock);
                    }

                    // Recompute distribution of predecessor state for further refinement.
                    recomputeChoiceDistributionsOfState(predecessorState);
                }
            }
        });
    }
}

template<typename ModelType>
storm::storage::Distribution<typename ModelType::ValueType, storm::storage::sparse::state_type>
EpsilonStableAbstractionDecomposition<ModelType>::computeEnhancedDistribution(
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

        auto lower = this->comparator.isLess(firstInterval.lower(), secondInterval.lower()) ? firstInterval.lower() : secondInterval.lower();
        auto upper = this->comparator.isLess(secondInterval.upper(), firstInterval.upper()) ? firstInterval.upper() : secondInterval.upper();

        enhancedDistribution.addProbability(key, ValueType(lower, carl::BoundType::WEAK, upper, carl::BoundType::WEAK));
    }

    return storm::utility::interval::makeDistributionProbabilistic(enhancedDistribution, this->comparator);
}

template<typename ModelType>
storm::IntervalBaseType<typename ModelType::ValueType> EpsilonStableAbstractionDecomposition<ModelType>::computeDeltaForState(
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& stateDistribution,
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& enhancedDistribution) {
    storm::IntervalBaseType<ValueType> delta = storm::utility::zero<IntervalBaseType<ValueType>>();

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
bool EpsilonStableAbstractionDecomposition<ModelType>::shouldBuildQuotient() const {
    return this->options.buildQuotient;
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::buildQuotientFromPartition() {
    // In order to create the quotient model, we need to construct
    // (a) the new transition matrix,
    // (b) the new labeling,
    // (c) the new reward structures.

    // Prepare a matrix builder for (a).
    // storm::storage::SparseMatrixBuilder<ValueType> builder(this->partition.getNumberOfBlocks(), this->partition.getNumberOfBlocks());
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
            auto choiceRangeOfRepresentativeState = getChoiceRangeOfState(representativeState);
            auto numberOfChoicesInBlock = choiceRangeOfRepresentativeState.second - choiceRangeOfRepresentativeState.first;
            for (uint_fast64_t choiceOffset = 0; choiceOffset < numberOfChoicesInBlock; choiceOffset++) {
                // Start block distribution for choice over final partition with the first state of the block.
                auto blockStateIterator = block.begin();
                auto currentBlockState = *blockStateIterator;
                auto choiceRangeOfCurrentBlockState = getChoiceRangeOfState(currentBlockState);
                auto blockChoiceDistribution = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();

                // Compute distribution of 'currentBlockState' over the quotient partition with current choice.
                for (auto const& e : this->model.getTransitionMatrix().getRow(choiceRangeOfCurrentBlockState.first + choiceOffset)) {
                    blockChoiceDistribution.addProbability(blocksMapping.at(this->partition.getBlockOfElement(e.getColumn()).front()), e.getValue());
                }

                // Now we enhance the block distribution with the intervals of every other state in the block.
                ++blockStateIterator;  // Go to next state of the block.
                while (blockStateIterator != block.end()) {
                    currentBlockState = *blockStateIterator;
                    choiceRangeOfCurrentBlockState = getChoiceRangeOfState(currentBlockState);

                    auto distributionOfCurrentBlockState = storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>();
                    for (auto const& e : this->model.getTransitionMatrix().getRow(choiceRangeOfCurrentBlockState.first + choiceOffset)) {
                        distributionOfCurrentBlockState.addProbability(blocksMapping.at(this->partition.getBlockOfElement(e.getColumn()).front()),
                                                                       e.getValue());
                    }
                    blockChoiceDistribution = computeEnhancedDistribution(blockChoiceDistribution, distributionOfCurrentBlockState);
                    ++blockStateIterator;
                }

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

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::debugFindMergeableFinalBlocks(double epsilon) {
    std::vector<storm::storage::stateminimization::Partition::Block> blocks;
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
bool EpsilonStableAbstractionDecomposition<ModelType>::canMergeBlocksForDebug(std::span<uint64_t const> blockA, std::span<uint64_t const> blockB,
                                                                              double epsilon) {
    std::vector<storm::storage::sparse::state_type> states;
    states.insert(states.end(), blockA.begin(), blockA.end());
    states.insert(states.end(), blockB.begin(), blockB.end());

    if (states.empty()) {
        return true;
    }

    auto candidate = originalChoiceDistributions[states.front()];

    for (std::size_t i = 1; i < states.size(); ++i) {
        candidate = computeEnhancedDistribution(candidate, originalChoiceDistributions[states[i]]);
    }

    bool ok = true;
    for (auto s : states) {
        auto delta = computeDeltaForState(originalChoiceDistributions[s], candidate);

        // std::cout << "merge-debug state " << s << ", delta = " << delta << ", 2delta = " << (2 * delta) << ", epsilon = " << epsilon << std::endl;

        if (2 * delta > epsilon) {
            ok = false;
        }
    }

    return ok;
}

// TODO: Add methods like this in a STORM_ASSERT_LOG?
template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::validateEpsilonStability(double epsilon) {
    bool valid = true;

    this->partition.forEachBlock([&](auto const& block) {
        if (block.empty() || this->absorbingBlocks.contains(block.front())) {
            return;
        }

        auto candidate = originalChoiceDistributions[block.front()];

        for (auto it = std::next(block.begin()); it != block.end(); ++it) {
            candidate = computeEnhancedDistribution(candidate, originalChoiceDistributions[*it]);
        }

        for (auto state : block) {
            auto delta = computeDeltaForState(originalChoiceDistributions[state], candidate);

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
void EpsilonStableAbstractionDecomposition<ModelType>::splitInitialPartitionBasedOnLabelsAndRewards() {
    for (auto const& label : this->options.respectedAtomicPropositions.value()) {
        if (label == "init") {
            continue;
        }

        auto labelledStates = this->model.getStates(label);
        this->partition.forEachBlock([this, &labelledStates](auto const& block) {
            this->partition.splitBlockByPredicate(block, [&labelledStates](auto const& e) { return labelledStates.get(e); });
        });
    }

    // If the model has state rewards, we need to consider them, because otherwise reward properties are not
    // preserved.
    if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
        // TODO: Check if this is implemented correctly
        this->splitInitialPartitionBasedOnRewards();
    }
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::splitInitialPartitionBasedOnActionSets() {
    std::vector<uint_fast64_t> actionIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    this->partition.forEachBlock([this, &actionIndices](auto const& block) {
        this->partition.splitBlockByOrder(block, [&actionIndices](auto const& s, auto const& t) {
            // Split by number of enabled actions.
            // TODO: Do we store the action labels?
            auto numberOfActionsS = (actionIndices[s + 1] - actionIndices[s]);
            auto numberOfActionsT = (actionIndices[t + 1] - actionIndices[t]);
            return numberOfActionsS < numberOfActionsT;
        });
    });
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::splitInitialPartitionBasedOnRewards() {
    auto const& rewardModel = this->model.getUniqueRewardModel();
    if (rewardModel.hasStateRewards()) {
        this->splitInitialPartitionBasedOnRewards(rewardModel.getStateRewardVector());
    }
    if (rewardModel.hasStateActionRewards()) {
        if (this->model.isNondeterministicModel()) {
            std::vector<std::set<ValueType>> actionRewards;
            actionRewards.reserve(this->model.getNumberOfStates());
            for (storm::storage::sparse::state_type state = 0; state < this->model.getNumberOfStates(); ++state) {
                std::set<ValueType> rewardsAtState;
                for (auto choice = this->model.getTransitionMatrix().getRowGroupIndices()[state];
                     choice < this->model.getTransitionMatrix().getRowGroupIndices()[state + 1]; ++choice) {
                    rewardsAtState.insert(rewardModel.getStateActionReward(choice));
                }
                actionRewards.push_back(std::move(rewardsAtState));
            }
            this->splitInitialPartitionBasedOnActionRewards(actionRewards);
        } else {
            this->splitInitialPartitionBasedOnRewards(rewardModel.getStateActionRewardVector());
        }
    }
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::splitInitialPartitionBasedOnRewards(std::vector<ValueType> const& rewardVector) {
    this->partition.forEachBlock([this, &rewardVector](auto const& block) {
        this->partition.splitBlockByOrder(block, [&rewardVector](auto const& a, auto const& b) { return rewardVector[a] < rewardVector[b]; });
    });
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::splitInitialPartitionBasedOnActionRewards(std::vector<std::set<ValueType>> const& actionRewards) {
    this->partition.forEachBlock([this, &actionRewards](auto const& block) {
        this->partition.splitBlockByOrder(block, [&actionRewards](auto const& a, auto const& b) { return actionRewards[a] < actionRewards[b]; });
    });
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::initializeChoiceDistributions() {
    auto const& matrix = this->model.getTransitionMatrix();

    originalChoiceDistributions.clear();
    originalChoiceDistributions.resize(matrix.getRowCount());

    for (auto choice = 0; choice < matrix.getRowCount(); ++choice) {
        recomputeChoiceDistribution(choice);
    }
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::recomputeChoiceDistribution(uint_fast64_t choice) {
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> distribution;

    for (auto const& e : this->model.getTransitionMatrix().getRow(choice)) {
        distribution.addProbability(this->partition.getBlockOfElement(e.getColumn()).front(), e.getValue());
    }

    originalChoiceDistributions[choice] = storm::utility::interval::makeDistributionProbabilistic(std::move(distribution), this->comparator);
}

template<typename ModelType>
std::span<storm::storage::Distribution<typename ModelType::ValueType, storm::storage::sparse::state_type> const>
EpsilonStableAbstractionDecomposition<ModelType>::getChoiceDistributionsOfState(storm::storage::sparse::state_type state) const {
    auto const [firstChoice, lastChoice] = getChoiceRangeOfState(state);

    return std::span<storm::storage::Distribution<typename ModelType::ValueType, storm::storage::sparse::state_type> const>(
        originalChoiceDistributions.data() + firstChoice, lastChoice - firstChoice);
}

template<typename ModelType>
void EpsilonStableAbstractionDecomposition<ModelType>::recomputeChoiceDistributionsOfState(storm::storage::sparse::state_type state) {
    auto const [firstChoice, lastChoice] = getChoiceRangeOfState(state);

    for (auto choice = firstChoice; choice < lastChoice; ++choice) {
        recomputeChoiceDistribution(choice);
    }
}

template<typename ModelType>
std::pair<std::uint_fast64_t, std::uint_fast64_t> EpsilonStableAbstractionDecomposition<ModelType>::getChoiceRangeOfState(
    storm::storage::sparse::state_type state) const {
    auto const& rowGroupIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    return {rowGroupIndices[state], rowGroupIndices[state + 1]};
}

template class EpsilonStableAbstractionDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;
template class EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>>;

template class EpsilonStableAbstractionDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>;
template class EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::RationalInterval>>;
}  // namespace abstraction
}  // namespace storage
}  // namespace storm
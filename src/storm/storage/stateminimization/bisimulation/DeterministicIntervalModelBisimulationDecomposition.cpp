#include "DeterministicIntervalModelBisimulationDecomposition.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::DeterministicIntervalModelBisimulationDecomposition(
    const ModelType& model, const typename BisimulationDecomposition<ModelType>::BisimulationOptions& options)
    : storm::storage::BisimulationDecomposition<ModelType>(model, options),
      probabilitiesToCurrentSplitter(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      probabilitiesToOtherBlocks(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      touchedProbabilitiesToSplitter(model.getNumberOfStates(), false) {}

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
    storm::storage::stateminimization::Partition::Block splitterBlock, std::deque<typename stateminimization::Partition::Block>& splitterQueue,
    stateminimization::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::stateminimization::Partition::BlockSet blocksToSplit;

    for (auto currentState : splitterBlock) {
        // TODO: Compute the transition intervals for each choice/action.
        for (const auto& predecessorEntry : this->backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);
            auto intervalTransitionProbability = predecessorEntry.getValue();

            if (!possiblyNeedsRefinement(predecessorBlock)) {
                continue;
            }

            // Compute interval going to the splitter, i.e., I(s, C), and initialize interval leading to all other blocks, i.e., I(s, \Pi \setminus C).
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
        // TODO: Compute the transition intervals for each choice/action.
        for (const auto& entry : this->model.getTransitionMatrix().getRow(predecessorState)) {
            auto targetState = entry.getColumn();

            if (!this->partition.contains(splitterBlock, targetState)) {
                probabilitiesToOtherBlocks[predecessorState] += entry.getValue();
            }
        }
    }

    for (auto predecessorBlockToSplit : blocksToSplit) {
        // First split the block by whether it is a predecessor of the splitter block or not.
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
                // TODO: Compare states for each possible choice/action.

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
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::buildQuotientFromPartition() {
    // TODO: Generalize the quotient construction to nondeterministic models (cf. EpsilonStableAbstractionDecomposition.cpp).
    // checkCurrentPartitionByExactFeasibleIntervals();

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

// TODO: Check this in an assert or during tests?
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

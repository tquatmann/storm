#include "storm/storage/stateminimization/bisimulation/DeterministicModelBisimulationDecomposition.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
DeterministicModelBisimulationDecomposition<ModelType>::DeterministicModelBisimulationDecomposition(
    ModelType const& model, typename BisimulationDecomposition<ModelType>::BisimulationOptions const& options)
    : BisimulationDecomposition<ModelType>(model, options),
      probabilitiesToCurrentSplitter(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      probabilitiesToOtherBlocks(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      touchedProbabilitiesToSplitter(model.getNumberOfStates(), false) {
    // Intentionally left empty.
}

template<typename ModelType>
std::pair<storm::storage::BitVector, storm::storage::BitVector> DeterministicModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    return storm::utility::graph::performProb01(this->backwardTransitions, this->options.phiStates.value(), this->options.psiStates.value());
}

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::splitOffDivergentStates() {
    // TODO: implement
}

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::initializeSilentProbabilities() {
    // TODO: implement
}

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::initializeWeakDtmcBisimulation() {
    // If we are creating the initial partition for weak bisimulation on DTMCs, we need to (a) split off all
    // divergent states of each initial block and (b) initialize the vector of silent probabilities.
    this->splitOffDivergentStates();
    this->initializeSilentProbabilities();
}

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::postProcessInitialPartition() {
    // TODO: implement
}

template<typename ModelType>
typename DeterministicModelBisimulationDecomposition<ModelType>::ValueType const&
DeterministicModelBisimulationDecomposition<ModelType>::getProbabilityToSplitter(storm::storage::sparse::state_type const& state) const {
    return probabilitiesToCurrentSplitter[state];
}

template<typename ModelType>
bool DeterministicModelBisimulationDecomposition<ModelType>::isSilent(storm::storage::sparse::state_type const& state) const {
    return this->comparator.isOne(silentProbabilities[state]);
}

template<typename ModelType>
bool DeterministicModelBisimulationDecomposition<ModelType>::hasNonZeroSilentProbability(storm::storage::sparse::state_type const& state) const {
    return !this->comparator.isZero(silentProbabilities[state]);
}

template<typename ModelType>
typename DeterministicModelBisimulationDecomposition<ModelType>::ValueType DeterministicModelBisimulationDecomposition<ModelType>::getSilentProbability(
    storm::storage::sparse::state_type const& state) const {
    return silentProbabilities[state];
}

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(
    std::span<uint64_t const> splitterBlock, std::deque<typename stateminimization::Partition::Block>& splitterQueue,
    stateminimization::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::stateminimization::Partition::BlockSet blocksToSplit;
    // std::cout << "Performing standard bisimulation!" << std::endl;

    // std::fill(probabilitiesToCurrentSplitter.begin(), probabilitiesToCurrentSplitter.end(), storm::utility::zero<ValueType>());
    for (auto currentState : splitterBlock) {
        // Compute probability to enter splitter block for each predecessor
        for (const auto& predecessorEntry : this->backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);
            auto transitionProbability = predecessorEntry.getValue();

            if (!possiblyNeedsRefinement(predecessorBlock)) {
                continue;
            }

            if (touchedProbabilitiesToSplitter.get(predecessorState)) {
                probabilitiesToCurrentSplitter[predecessorState] += transitionProbability;
            } else {
                probabilitiesToCurrentSplitter[predecessorState] = transitionProbability;
                touchedProbabilitiesToSplitter.set(predecessorState, true);
            }

            // Remember which blocks contain predecessors to split them w.r.t. the splitter afterwards
            blocksToSplit.insert(predecessorBlock);
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
                return this->comparator.isLess(probabilitiesToCurrentSplitter[a], probabilitiesToCurrentSplitter[b]);
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

// template<typename ModelType>
// void DeterministicModelBisimulationDecomposition<ModelType>::refinePredecessorBlockOfSplitterStrong()
//     // TODO: implement
// }
//
// template<typename ModelType>
// void DeterministicModelBisimulationDecomposition<ModelType>::refinePredecessorBlocksOfSplitterStrong()
//     // TODO: implement
// }

template<typename ModelType>
bool DeterministicModelBisimulationDecomposition<ModelType>::possiblyNeedsRefinement(std::span<uint64_t const> block) const {
    return block.size() > 1 && !this->absorbingBlocks.contains(block.front());
}

//
// template<typename ModelType>
// void DeterministicModelBisimulationDecomposition<ModelType>::increaseProbabilityToSplitter(storm::storage::sparse::state_type predecessor,
//                                                                                            bisimulation::Block<BlockDataType> const& predecessorBlock,
//                                                                                            ValueType const& value) {
//     STORM_LOG_TRACE("Increasing probability of " << predecessor << " to splitter by " << value << ".");
//     storm::storage::sparse::state_type predecessorPosition = this->partition.getPosition(predecessor);
//
//     // If the position of the state is to the right of marker1, we have not seen it before.
//     if (predecessorPosition >= predecessorBlock.data().marker1()) {
//         // Then, we just set the value.
//         probabilitiesToCurrentSplitter[predecessor] = value;
//     } else {
//         // If the state was seen as a predecessor before, we add the value to the existing value.
//         probabilitiesToCurrentSplitter[predecessor] += value;
//     }
// }

template<typename ModelType>
void DeterministicModelBisimulationDecomposition<ModelType>::buildQuotientFromPartition() {
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

            // Add all of the selected atomic propositions that hold in the representative state to the state
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
                    probIterator->second += getTransitionValue(entry, representativeState);
                } else {
                    blockProbability[targetBlock] = getTransitionValue(entry, representativeState);
                    ;
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
    this->quotient = std::make_shared<ModelType>(builder.build(), std::move(newLabeling), std::move(rewardModels));
}

template<typename ModelType>
DeterministicModelBisimulationDecomposition<ModelType>::ValueType DeterministicModelBisimulationDecomposition<ModelType>::getTransitionValue(
    storm::storage::MatrixEntry<storm::storage::sparse::state_type, ValueType> const& matrixEntry, storm::storage::sparse::state_type state) const {
    if constexpr (std::is_same_v<ModelType, storm::models::sparse::Ctmc<typename ModelType::ValueType>>) {
        auto transitionValue = matrixEntry.getValue();
        // TODO: enable when removing CTMC rate matrix
        // transitionValue *= this->model.getExitRateVector().at(state);
        return transitionValue;
    } else {
        STORM_LOG_ASSERT(this->model.isDiscreteTimeModel(), "Unhandled model type");
        return matrixEntry.getValue();
    }
}

template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Dtmc<double>>;
template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Ctmc<double>>;

template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalNumber>>;
template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalNumber>>;

template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalFunction>>;
template class DeterministicModelBisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalFunction>>;
}  // namespace storage
}  // namespace storm

#include "storm/storage/stateminimization/bisimulation/BisimulationDecomposition.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationOptions::BisimulationOptions(ModelType const& model, storm::logic::Formula const& formula)
    : BaseDecomposition<ModelType>::BaseOptions(model, formula) {
    // Intentionally left empty.
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationOptions::BisimulationOptions(ModelType const& model,
                                                                               std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas)
    : BaseDecomposition<ModelType>::BaseOptions(model, formulas) {
    // Intentionally left empty.
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationOptions::BisimulationOptions() : BaseDecomposition<ModelType>::BaseOptions(), type(BisimulationType::Strong) {
    // Intentionally left empty.
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationDecomposition(ModelType const& model, BisimulationOptions const& options)
    : BisimulationDecomposition(model, model.getBackwardTransitions(), options) {
    // Intentionally left empty.
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationDecomposition(ModelType const& model, storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                                BisimulationOptions const& options)
    : BaseDecomposition<ModelType>(model, backwardTransitions), options(options), comparator(options.getTolerance()) {
    STORM_LOG_THROW(!options.getKeepRewards() || !model.hasRewardModel() || model.hasUniqueRewardModel(), storm::exceptions::IllegalFunctionCallException,
                    "Bisimulation currently only supports models with at most one reward model.");
    STORM_LOG_THROW(!options.getKeepRewards() || !model.hasRewardModel() || !model.getUniqueRewardModel().hasTransitionRewards(),
                    storm::exceptions::IllegalFunctionCallException,
                    "Bisimulation is currently supported for models with state or action rewards only. Consider converting the transition rewards to state "
                    "rewards (via suitable function calls).");
    STORM_LOG_THROW(options.getType() != BisimulationType::Weak || !options.isBounded(), storm::exceptions::IllegalFunctionCallException,
                    "Weak bisimulation cannot preserve bounded properties.");
    STORM_LOG_THROW(options.getType() != BisimulationType::Weak || !options.isDiscounted(), storm::exceptions::IllegalFunctionCallException,
                    "Weak bisimulation cannot preserve discounted properties.");

    // Fix the respected atomic propositions if they were not explicitly given.
    if (!this->options.respectedAtomicPropositions) {
        this->options.respectedAtomicPropositions = model.getStateLabeling().getLabels();
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::refinePartition() {
    this->performPartitionRefinement();
    // this->performSignatureRefinement();
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::computeInitialPartition() {
    if (false /*options.measureDrivenInitialPartition*/) {
        std::cout << "Initializing measure driven!" << std::endl;
        STORM_LOG_THROW(this->options.phiStates, storm::exceptions::InvalidOptionException,
                        "Unable to compute measure-driven initial partition without phi states.");
        STORM_LOG_THROW(this->options.psiStates, storm::exceptions::InvalidOptionException,
                        "Unable to compute measure-driven initial partition without psi states.");
        this->initializeMeasureDrivenPartition();
    } else {
        this->initializeLabelBasedPartition();
    }
}

template<typename ModelType>
bool BisimulationDecomposition<ModelType>::shouldBuildQuotient() const {
    return this->options.buildQuotient;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performPartitionRefinement() {
    std::deque<typename stateminimization::Partition::Block> splitterQueue;
    stateminimization::Partition::BlockSet enqueuedSplitterBlocks;

    // Initially, add all current blocks to the queue.
    this->partition.forEachBlock([&](auto const& block) {
        splitterQueue.push_back(block);
        enqueuedSplitterBlocks.insert(block);
    });

    // Then perform the actual splitting until there are no more splitters.
    uint_fast64_t iterations = 0;
    while (!splitterQueue.empty()) {
        ++iterations;

        auto splitterBlock = splitterQueue.front();
        splitterQueue.pop_front();
        enqueuedSplitterBlocks.erase(splitterBlock);

        refinePartitionBasedOnSplitter(splitterBlock, splitterQueue, enqueuedSplitterBlocks);
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performSignatureRefinement() {
    // Insert all blocks into the queue for refinement
    std::deque<typename stateminimization::Partition::Block> blocksQueue;

    storm::storage::stateminimization::Partition::BlockSet enqueuedSplitterBlocks;
    // Initially, add all current blocks to the queue
    this->partition.forEachBlock([&](auto const& block) {
        // TODO: maybe one has to handle the absorbing blocks differently for weak bisimulation, i.e., that they are still enqueued here and handled
        // differently while splitting based on the computed signatures
        if (!this->absorbingBlocks.contains(block.front())) {
            blocksQueue.push_back(block);
            enqueuedSplitterBlocks.insert(block);
        }
    });

    std::vector<size_t> stateToSignature(this->backwardTransitions.getColumnCount(), 0);
    std::vector<storm::storage::sparse::state_type> statesWithInvalidSignature;

    // Refine the partition as long as the queue is not empty
    uint_fast64_t iterations = 0;
    uint_fast64_t noSplitCounter = 0;
    while (!blocksQueue.empty()) {
        ++iterations;

        auto blockToRefine = blocksQueue.back();
        blocksQueue.pop_back();
        enqueuedSplitterBlocks.erase(blockToRefine);

        // Detect if splitting is necessary
        size_t firstSignature = 0;
        bool hasMultipleSignatures = false;
        bool isFirstState = true;

        // Map states to their signature
        for (auto state : blockToRefine) {
            // Only compute the state signature if it was invalidated
            if (stateToSignature[state] == 0) {
                auto signatureHash = computeStateSignatureHash(state);
                stateToSignature[state] = signatureHash;
            }

            // Compare signatures
            if (isFirstState) {
                firstSignature = stateToSignature[state];
                isFirstState = false;
            } else if (stateToSignature[state] != firstSignature) {
                hasMultipleSignatures = true;
            }
        }

        // Skip splitting if all states have the same signature
        if (!hasMultipleSignatures) {
            noSplitCounter++;
            continue;
        }

        auto oldFront = blockToRefine.front();
        // split blocks according to their state signatures, if possible
        auto wasSplit = this->partition.splitBlockByOrder(
            blockToRefine, [&stateToSignature](auto const& a, auto const& b) { return stateToSignature.at(a) < stateToSignature.at(b); });

        if (wasSplit) {
            auto newFront = blockToRefine.front();
            if (newFront != oldFront) {
                std::cout << "[LOG] Block front changed: was " << oldFront << ", now " << newFront << std::endl;
            }
            this->partition.forEachSubBlock(blockToRefine, [this, &blocksQueue, &statesWithInvalidSignature, &enqueuedSplitterBlocks](auto const& block) {
                // TODO: If representative state is already on queue, then don't add it
                if (block.size() > 1 && !enqueuedSplitterBlocks.contains(block)) {
                    blocksQueue.push_back(block);
                    enqueuedSplitterBlocks.insert(block);
                }

                for (auto state : block) {
                    for (auto& transition : this->backwardTransitions.getRow(state)) {
                        auto predecessorState = transition.getColumn();
                        auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);

                        // place target block on queue only if it is not already there
                        if (predecessorBlock.size() > 1 && !enqueuedSplitterBlocks.contains(predecessorBlock)) {
                            blocksQueue.push_back(predecessorBlock);
                            enqueuedSplitterBlocks.insert(predecessorBlock);
                        }

                        // remember which states have an invalid signature now
                        statesWithInvalidSignature.emplace_back(predecessorState);
                    }
                }
            });
        }

        // invalidate signatures of affected states
        for (auto currentState : statesWithInvalidSignature) {
            stateToSignature[currentState] = 0;
        }
        statesWithInvalidSignature.clear();

        if (!wasSplit) {
            noSplitCounter++;
        }

        if (storm::utility::resources::isTerminate()) {
            std::cout << "Performed " << iterations << " iterations of partition refinement before abort.\n";
            STORM_LOG_THROW(false, storm::exceptions::AbortException, "Aborted in bisimulation computation.");
            break;
        }
    }

    std::cout << "Finished refinement after " << iterations << " iterations." << std::endl;
    std::cout << "Attempt to split block failed " << noSplitCounter << " times." << std::endl;
}

template<typename ModelType>
std::size_t BisimulationDecomposition<ModelType>::computeStateSignatureHash(storm::storage::sparse::state_type state) const {
    // TODO: Interestingly, the boost hash function seems to create collisions, hence the quotient is not minimal per default.
    // TODO: When using --exact, the quotient gets calculated correctly. Thus, Boost does not guarantee a sufficient
    // TODO: precision for doubles
    return computeStateSignature(state, this->partition).computeHash();
    // TODO: Using --exact with the string based hash representation does not lead to a correct quotient -> investigate
    // return std::hash<std::string>{}(computeStateSignature(state, partition).toString());
}

template<typename ModelType>
storm::storage::bisimulation::Signature<typename ModelType::ValueType> BisimulationDecomposition<ModelType>::computeStateSignature(
    storm::storage::sparse::state_type state, storm::storage::stateminimization::Partition const& currentPartition) const {
    storm::storage::bisimulation::Signature<typename ModelType::ValueType> signature;

    for (auto entry : this->model.getTransitionMatrix().getRow(state)) {
        // std::cout << "Prob for state " << state << " to reach target state " << entry.getColumn() << ": " << entry.getValue() << std::endl;
        auto targetBlock = this->partition.getBlockOfElement(entry.getColumn());  // column marks the id of the target state
        signature.addBlockProbability(targetBlock.front(), entry.getValue());
    }

    return signature;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnRewards() {
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
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnRewards(std::vector<ValueType> const& rewardVector) {
    this->partition.forEachBlock([this, &rewardVector](auto const& block) {
        this->partition.splitBlockByOrder(block, [&rewardVector](auto const& a, auto const& b) { return rewardVector[a] < rewardVector[b]; });
    });
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnActionRewards(std::vector<std::set<ValueType>> const& actionRewards) {
    this->partition.forEachBlock([this, &actionRewards](auto const& block) {
        this->partition.splitBlockByOrder(block, [&actionRewards](auto const& a, auto const& b) { return actionRewards[a] < actionRewards[b]; });
    });
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::initializeLabelBasedPartition() {
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

    std::cout << "Number of blocks after initial partitioning by labels: " << this->partition.getNumberOfBlocks() << std::endl;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::initializeMeasureDrivenPartition() {
    std::pair<storm::storage::BitVector, storm::storage::BitVector> statesWithProbability01 = this->getStatesWithProbability01();

    std::optional<storm::storage::sparse::state_type> representativePsiState;
    if (!this->options.psiStates.value().empty()) {
        representativePsiState = *this->options.psiStates.value().begin();
    }

    auto prob0States = statesWithProbability01.first;
    auto prob1States = statesWithProbability01.second;

    // TODO: These blocks should be marked as absorbing
    auto initialBlock = this->partition.getBlockOfElement(0);
    if (!prob0States.empty()) {
        auto [notProb0Block, prob0Block] =
            this->partition.splitBlockByPredicate(initialBlock, [&prob0States](auto const& state) { return prob0States.get(state); });
        this->absorbingBlocks.emplace(prob0Block.front(), prob0Block.front());

        if (!prob1States.empty()) {
            auto [notProb1Block, prob1Block] =
                this->partition.splitBlockByPredicate(notProb0Block, [&prob1States](auto const& state) { return prob1States.get(state); });
            this->absorbingBlocks.emplace(prob1Block.front(), representativePsiState.value());
        }
    } else if (!prob1States.empty()) {
        auto [notProb1Block, prob1Block] =
            this->partition.splitBlockByPredicate(initialBlock, [&prob1States](auto const& state) { return prob1States.get(state); });
        this->absorbingBlocks.emplace(prob1Block.front(), representativePsiState.value());

        if (!prob0States.empty()) {
            auto [notProb0Block, prob0Block] =
                this->partition.splitBlockByPredicate(notProb1Block, [&prob0States](auto const& state) { return prob0States.get(state); });
            this->absorbingBlocks.emplace(prob0Block.front(), prob0Block.front());
        }
    }

    // If the model has state rewards, we need to consider them, because otherwise reward properties are not
    // preserved.
    if (this->options.getKeepRewards() && this->model.hasRewardModel()) {
        this->splitInitialPartitionBasedOnRewards();
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::initialize() {
    // Intentionally left empty.
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::extractDecompositionBlocks() {
    // TODO: implement
}

template class BisimulationDecomposition<storm::models::sparse::Dtmc<double>>;
template class BisimulationDecomposition<storm::models::sparse::Ctmc<double>>;
template class BisimulationDecomposition<storm::models::sparse::Mdp<double>>;

template class BisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalNumber>>;
template class BisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalNumber>>;
template class BisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalNumber>>;

template class BisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalFunction>>;
template class BisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalFunction>>;
template class BisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalFunction>>;

template class storm::storage::BisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;
template class storm::storage::BisimulationDecomposition<storm::models::sparse::Ctmc<storm::Interval>>;
template class storm::storage::BisimulationDecomposition<storm::models::sparse::Mdp<storm::Interval>>;

template class storm::storage::BisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>;
template class storm::storage::BisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalInterval>>;
template class storm::storage::BisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalInterval>>;
}  // namespace storage
}  // namespace storm

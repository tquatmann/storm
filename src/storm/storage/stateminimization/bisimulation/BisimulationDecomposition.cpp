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
    : BaseDecomposition<ModelType>(model, backwardTransitions, options.getTolerance()), options(options) {
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
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::computeInitialPartition() {
    // TODO: re-enable measureDriven
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

    if (this->model.isNondeterministicModel()) {
        this->splitInitialPartitionBasedOnActionSets();
    }

    std::cout << "Size after initial partition: " << this->partition.getNumberOfBlocks() << std::endl;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnActionSets() {
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
bool BisimulationDecomposition<ModelType>::shouldBuildQuotient() const {
    return this->options.buildQuotient;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performPartitionRefinement() {
    std::deque<typename stateminimization::Partition::Block> splitterQueue; // todo: not used
    stateminimization::Partition::BlockSet enqueuedSplitterBlocks;

    // Initially, add all current blocks to the queue.
    this->partition.forEachBlock([&](auto const& block) {
        enqueuedSplitterBlocks.insert(block);
    });

    // Then perform the actual splitting until there are no more splitters.
    uint64_t iterations = 0;
    while (!enqueuedSplitterBlocks.empty()) {
        ++iterations;

        auto splitterBlock = *enqueuedSplitterBlocks.begin();
        // std::cout << "Iteration " << iterations << ", number of blocks in queue: " << enqueuedSplitterBlocks.size() << " with sizes ";
        //  for (auto const& block : enqueuedSplitterBlocks) {
        //      std::cout << block.size() << " ";
        //  }
        // std::cout << "Picked block with " << splitterBlock.size() << " states." << std::endl;
        enqueuedSplitterBlocks.erase(enqueuedSplitterBlocks.begin());
        if (this->partition.isProperSuperBlock(splitterBlock)) {
            // If the splitter block is a proper super block, we can skip it as we already have checked all sub-blocks TODO: apparently this is not true, why?
            // continue;
        }

        refinePartitionBasedOnSplitter(splitterBlock, splitterQueue, enqueuedSplitterBlocks);
    }
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

    if (this->options.isActionSensitive()) {
        this->splitInitialPartitionBasedOnActionSets();
    }
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

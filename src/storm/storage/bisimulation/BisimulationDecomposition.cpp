#include "storm/storage/bisimulation/BisimulationDecomposition.h"

#include <chrono>

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/exceptions/AbortException.h"
#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/exceptions/InvalidOptionException.h"

#include "storm/logic/FormulaInformation.h"
#include "storm/logic/FragmentSpecification.h"

#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/modelchecker/propositional/SparsePropositionalModelChecker.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/CoreSettings.h"

#include "Signature.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/macros.h"

#include <optional>  // add this
#include "BucketKey.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
BisimulationDecomposition<ModelType>::Options::Options(ModelType const& model, storm::logic::Formula const& formula) : Options() {
    this->preserveSingleFormula(model, formula);
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::Options::Options(ModelType const& model,
                                                                      std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas)
    : Options() {
    if (formulas.empty()) {
        this->respectedAtomicPropositions = model.getStateLabeling().getLabels();
        this->keepRewards = true;
    }
    if (formulas.size() == 1) {
        this->preserveSingleFormula(model, *formulas.front());
    } else {
        for (auto const& formula : formulas) {
            preserveFormula(*formula);
        }
    }
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::Options::Options()
    : measureDrivenInitialPartition(false),
      phiStates(),
      psiStates(),
      respectedAtomicPropositions(),
      buildQuotient(true),
      keepRewards(false),
      type(BisimulationType::Strong),
      bounded(false),
      usesEpsilon(false) {
    // Intentionally left empty.
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::Options::preserveFormula(storm::logic::Formula const& formula) {
    // Disable the measure driven initial partition.
    measureDrivenInitialPartition = false;
    phiStates = boost::none;
    psiStates = boost::none;

    // Retrieve information about formula.
    storm::logic::FormulaInformation info = formula.info();

    // Preserve rewards if necessary.
    keepRewards = keepRewards || info.containsRewardOperator() || info.containsRewardBoundedFormula();

    // Preserve bounded properties if necessary.
    bounded = bounded || (info.containsBoundedUntilFormula() || info.containsNextFormula() || info.containsCumulativeRewardFormula());

    // Compute the relevant labels and expressions.
    this->addToRespectedAtomicPropositions(formula.getAtomicExpressionFormulas(), formula.getAtomicLabelFormulas());
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::Options::preserveSingleFormula(ModelType const& model, storm::logic::Formula const& formula) {
    // Retrieve information about formula.
    storm::logic::FormulaInformation info = formula.info();

    keepRewards = info.containsRewardOperator() || info.containsRewardBoundedFormula();

    // We need to preserve bounded properties iff the formula contains a bounded until or a next subformula.
    bounded = info.containsBoundedUntilFormula() || info.containsNextFormula() || info.containsCumulativeRewardFormula();

    // Compute the relevant labels and expressions.
    this->addToRespectedAtomicPropositions(formula.getAtomicExpressionFormulas(), formula.getAtomicLabelFormulas());

    // Check whether measure driven initial partition is possible and, if so, set it.
    this->checkAndSetMeasureDrivenInitialPartition(model, formula);
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::Options::checkAndSetMeasureDrivenInitialPartition(ModelType const& model,
                                                                                                            storm::logic::Formula const& formula) {
    std::shared_ptr<storm::logic::Formula const> newFormula = formula.asSharedPointer();

    if (formula.isProbabilityOperatorFormula()) {
        if (formula.asProbabilityOperatorFormula().hasOptimalityType()) {
            optimalityType = formula.asProbabilityOperatorFormula().getOptimalityType();
        } else if (formula.asProbabilityOperatorFormula().hasBound()) {
            storm::logic::ComparisonType comparisonType = formula.asProbabilityOperatorFormula().getComparisonType();
            if (comparisonType == storm::logic::ComparisonType::Less || comparisonType == storm::logic::ComparisonType::LessEqual) {
                optimalityType = OptimizationDirection::Maximize;
            } else {
                optimalityType = OptimizationDirection::Minimize;
            }
        }
        newFormula = formula.asProbabilityOperatorFormula().getSubformula().asSharedPointer();
    } else if (formula.isRewardOperatorFormula()) {
        if (formula.asRewardOperatorFormula().hasOptimalityType()) {
            optimalityType = formula.asRewardOperatorFormula().getOptimalityType();
        } else if (formula.asRewardOperatorFormula().hasBound()) {
            storm::logic::ComparisonType comparisonType = formula.asRewardOperatorFormula().getComparisonType();
            if (comparisonType == storm::logic::ComparisonType::Less || comparisonType == storm::logic::ComparisonType::LessEqual) {
                optimalityType = OptimizationDirection::Maximize;
            } else {
                optimalityType = OptimizationDirection::Minimize;
            }
        }
        newFormula = formula.asRewardOperatorFormula().getSubformula().asSharedPointer();
    }

    std::shared_ptr<storm::logic::Formula const> leftSubformula = std::make_shared<storm::logic::BooleanLiteralFormula>(true);
    std::shared_ptr<storm::logic::Formula const> rightSubformula;
    if (newFormula->isUntilFormula()) {
        leftSubformula = newFormula->asUntilFormula().getLeftSubformula().asSharedPointer();
        rightSubformula = newFormula->asUntilFormula().getRightSubformula().asSharedPointer();
        if (leftSubformula->isInFragment(storm::logic::propositional()) && rightSubformula->isInFragment(storm::logic::propositional())) {
            measureDrivenInitialPartition = true;
        }
    } else if (newFormula->isEventuallyFormula()) {
        rightSubformula = newFormula->asEventuallyFormula().getSubformula().asSharedPointer();
        if (rightSubformula->isInFragment(storm::logic::propositional())) {
            measureDrivenInitialPartition = true;
        }
    }

    if (measureDrivenInitialPartition) {
        storm::modelchecker::SparsePropositionalModelChecker<ModelType> checker(model);
        std::unique_ptr<storm::modelchecker::CheckResult> phiStatesCheckResult = checker.check(*leftSubformula);
        std::unique_ptr<storm::modelchecker::CheckResult> psiStatesCheckResult = checker.check(*rightSubformula);
        phiStates = phiStatesCheckResult->asExplicitQualitativeCheckResult().getTruthValuesVector();
        psiStates = psiStatesCheckResult->asExplicitQualitativeCheckResult().getTruthValuesVector();
    } else {
        optimalityType = boost::none;
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::Options::addToRespectedAtomicPropositions(
    std::vector<std::shared_ptr<storm::logic::AtomicExpressionFormula const>> const& expressions,
    std::vector<std::shared_ptr<storm::logic::AtomicLabelFormula const>> const& labels) {
    std::set<std::string> labelsToRespect;
    for (auto const& labelFormula : labels) {
        labelsToRespect.insert(labelFormula->getLabel());
    }
    for (auto const& expressionFormula : expressions) {
        labelsToRespect.insert(expressionFormula->toString());
    }
    if (!respectedAtomicPropositions) {
        respectedAtomicPropositions = labelsToRespect;
    } else {
        respectedAtomicPropositions.get().insert(labelsToRespect.begin(), labelsToRespect.end());
    }
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationDecomposition(ModelType const& model, Options const& options)
    : BisimulationDecomposition(model, model.getBackwardTransitions(), options) {
    // Intentionally left empty.
}

template<typename ModelType>
BisimulationDecomposition<ModelType>::BisimulationDecomposition(ModelType const& model,
                                                                               storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                                               Options const& options)
    : model(model), backwardTransitions(backwardTransitions), options(options), partition(model.getNumberOfStates()), comparator(), quotient(nullptr), absorbingBlocks() {
    STORM_LOG_THROW(!options.getKeepRewards() || !model.hasRewardModel() || model.hasUniqueRewardModel(), storm::exceptions::IllegalFunctionCallException,
                    "Bisimulation currently only supports models with at most one reward model.");
    STORM_LOG_THROW(!options.getKeepRewards() || !model.hasRewardModel() || !model.getUniqueRewardModel().hasTransitionRewards(),
                    storm::exceptions::IllegalFunctionCallException,
                    "Bisimulation is currently supported for models with state or action rewards only. Consider converting the transition rewards to state "
                    "rewards (via suitable function calls).");
    STORM_LOG_THROW(options.getType() != BisimulationType::Weak || !options.getBounded(), storm::exceptions::IllegalFunctionCallException,
                    "Weak bisimulation cannot preserve bounded properties.");

    // Fix the respected atomic propositions if they were not explicitly given.
    if (!this->options.respectedAtomicPropositions) {
        this->options.respectedAtomicPropositions = model.getStateLabeling().getLabels();
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::computeBisimulationDecomposition() {
    std::chrono::high_resolution_clock::time_point totalStart = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point initialPartitionStart = std::chrono::high_resolution_clock::now();
    // initialize the initial partition.
    if (false/*options.measureDrivenInitialPartition*/) {
        std::cout << "Initializing measure driven!" << std::endl;
        STORM_LOG_THROW(options.phiStates, storm::exceptions::InvalidOptionException, "Unable to compute measure-driven initial partition without phi states.");
        STORM_LOG_THROW(options.psiStates, storm::exceptions::InvalidOptionException, "Unable to compute measure-driven initial partition without psi states.");
        this->initializeMeasureDrivenPartition();
    } else {
        this->initializeLabelBasedPartition();
    }
    STORM_LOG_WARN_COND(partition.getNumberOfBlocks() > 1, "Initial partition consists only of a single block.");
    std::chrono::high_resolution_clock::duration initialPartitionTime = std::chrono::high_resolution_clock::now() - initialPartitionStart;

    this->initialize();

    std::chrono::high_resolution_clock::time_point refinementStart = std::chrono::high_resolution_clock::now();

    if (this->options.getUsesEpsilon()) {
        auto epsilon = this->options.getEpsilon();
        this->performEpsilonSignatureRefinementUsingCompleteLinkage(epsilon);
    } else {
        this->performPartitionRefinement();
        // this->performSignatureRefinement();
    }

    std::chrono::high_resolution_clock::duration refinementTime = std::chrono::high_resolution_clock::now() - refinementStart;

    std::chrono::milliseconds refinementTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(refinementTime);
    std::cout << "    * time for partitioning: " << refinementTimeInMilliseconds.count() << "ms\n";

    std::chrono::high_resolution_clock::time_point extractionStart = std::chrono::high_resolution_clock::now();
    this->extractDecompositionBlocks();
    std::chrono::high_resolution_clock::duration extractionTime = std::chrono::high_resolution_clock::now() - extractionStart;

    std::chrono::high_resolution_clock::time_point quotientBuildStart = std::chrono::high_resolution_clock::now();
    if (options.buildQuotient) {
        this->buildQuotient();
    }
    std::chrono::high_resolution_clock::duration quotientBuildTime = std::chrono::high_resolution_clock::now() - quotientBuildStart;

    std::chrono::high_resolution_clock::duration totalTime = std::chrono::high_resolution_clock::now() - totalStart;

    if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isShowStatisticsSet()) {
        std::chrono::milliseconds initialPartitionTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(initialPartitionTime);
        std::chrono::milliseconds refinementTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(refinementTime);
        std::chrono::milliseconds extractionTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(extractionTime);
        std::chrono::milliseconds quotientBuildTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(quotientBuildTime);
        std::chrono::milliseconds totalTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(totalTime);
        std::cout << "\nTime breakdown:\n";
        std::cout << "    * time for initial partition: " << initialPartitionTimeInMilliseconds.count() << "ms\n";
        std::cout << "    * time for partitioning: " << refinementTimeInMilliseconds.count() << "ms\n";
        std::cout << "    * time for extraction: " << extractionTimeInMilliseconds.count() << "ms\n";
        std::cout << "    * time for building quotient: " << quotientBuildTimeInMilliseconds.count() << "ms\n";
        std::cout << "------------------------------------------\n";
        std::cout << "    * total time: " << totalTimeInMilliseconds.count() << "ms\n\n";
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performPartitionRefinement() {
    std::deque<typename Partition::Block> splitterQueue;
    bisimulation::Partition::BlockSet enqueuedSplitterBlocks;

    // Initially, add all current blocks to the queue
    this->partition.forEachBlock([&](auto const& block) {
        splitterQueue.push_back(block);
        enqueuedSplitterBlocks.insert(block);
    });

    uint_fast64_t iterations = 0;
    while (!splitterQueue.empty()) {
        ++iterations;

        auto splitterBlock = splitterQueue.front();
        splitterQueue.pop_front();
        enqueuedSplitterBlocks.erase(splitterBlock);

        refinePartitionBasedOnSplitter(splitterBlock, splitterQueue, enqueuedSplitterBlocks);
    }

    std::cout << "Finished refinement after " << iterations << " iterations." << std::endl;
    std::cout << "Size of final partition " << partition.getNumberOfBlocks() << "." << std::endl;
}

// uint_fast64_t iterations = 0;
//
// uint_fast64_t numberOfBlocks;
// do {
//     numberOfBlocks = partition.getNumberOfBlocks();
//
//     partition.forEachBlock([this, &iterations] (auto const& splitterBlock) {
//         ++iterations;
//
//         refinePartitionBasedOnSplitter(splitterBlock);
//     });
// } while (partition.getNumberOfBlocks() > numberOfBlocks);
//
// std::cout << "Finished refinement after " << iterations << " iterations." << std::endl;


template<typename ModelType>
void BisimulationDecomposition<ModelType>::performEpsilonSignatureRefinement(double epsilon) {
    std::cout << "Computing epsilon-bisimulation with epsilon: " << epsilon << std::endl;

    if constexpr (!storm::IsIntervalType<ValueType>) {
        STORM_LOG_THROW(true, storm::exceptions::IllegalFunctionCallException, "Cannot compute epsilon-bisimulation on non-interval model!");
    }

    // Insert all blocks into the queue for refinement
    std::deque<typename Partition::Block> blocksQueue;

    storm::storage::Partition::BlockSet enqueuedSplitterBlocks;
    // Initially, add all current blocks to the queue
    this->partition.forEachBlock([&](auto const& block) {
        // TODO: maybe one has to handle the absorbing blocks differently for weak bisimulation, i.e., that they are still enqueued here and handled differently
        // while splitting based on the computed signatures
        if (!this->absorbingBlocks.contains(block.front())) {
            blocksQueue.push_back(block);
            enqueuedSplitterBlocks.insert(block);
        }
    });

    std::vector<std::optional<storm::storage::bisimulation::BucketKey>> stateToBucketKey(backwardTransitions.getColumnCount());
    std::vector<storm::storage::sparse::state_type> statesWithInvalidBucketKey;

    // Refine the partition as long as the queue is not empty
    uint_fast64_t iterations = 0;
    uint_fast64_t noSplitCounter = 0;
    while (!blocksQueue.empty()) {
        ++iterations;

        auto blockToRefine = blocksQueue.back();
        blocksQueue.pop_back();
        enqueuedSplitterBlocks.erase(blockToRefine);

        auto projectionFn = [&](std::uint64_t state, std::span<const std::uint64_t> block) -> std::pair<double,double> {
            if constexpr (storm::IsIntervalType<ValueType>) {
                double lower = 0.0;
                double upper = 0.0;

                auto intervalToBlock = storm::utility::zero<ValueType>();

                for (auto const& entry : model.getTransitionMatrix().getRow(state)) {
                    if (!partition.contains(block, entry.getColumn())) continue;
                    intervalToBlock += entry.getValue();
                }

                // normalize to [0,1]
                lower = std::min(1.0, intervalToBlock.lower());
                upper = std::min(1.0, intervalToBlock.upper());
                return {lower, upper};
            } else {
                STORM_LOG_THROW(true, storm::exceptions::IllegalFunctionCallException, "Cannot compute epsilon-bisimulation on non-interval model!");
            }
        };

        storm::storage::bisimulation::BucketKeyBuilder keyBuilder(partition, epsilon, projectionFn);
        // TODO: Compute intervals to each block for every state in blockToRefine
        // TODO: Create BucketKey for every state based on intervals
        for (auto state : blockToRefine) {
            if (!stateToBucketKey[state].has_value()) {
                stateToBucketKey[state] = keyBuilder.buildForState(state);
            }
        }

        // TODO: Split states according to BucketKey; invalidate intervals for all predecessors (and their corresponding blocks)
        auto wasSplit = this->partition.splitBlockByOrder(blockToRefine, [&stateToBucketKey]
                                                          (auto const& a, auto const& b) { return stateToBucketKey.at(a) < stateToBucketKey.at(b); });

        if (wasSplit) {
            this->partition.forEachSubBlock(blockToRefine, [this, &blocksQueue, &statesWithInvalidBucketKey, &enqueuedSplitterBlocks](auto const& block) {
                // TODO: If representative state is already on queue, then don't add it
                if (block.size() > 1 && !enqueuedSplitterBlocks.contains(block)) {
                    blocksQueue.push_back(block);
                    enqueuedSplitterBlocks.insert(block);
                }

                for (auto state : block) {
                    for (auto &transition: backwardTransitions.getRow(state)) {
                        auto predecessorState = transition.getColumn();
                        auto predecessorBlock = partition.getBlockOfElement(predecessorState);

                        // place target block on queue only if it is not already there
                        if (predecessorBlock.size() > 1 && !enqueuedSplitterBlocks.contains(predecessorBlock)) {
                            blocksQueue.push_back(predecessorBlock);
                            enqueuedSplitterBlocks.insert(predecessorBlock);
                        }

                        // remember which states have an invalid signature now
                        statesWithInvalidBucketKey.emplace_back(predecessorState);
                    }
                }
            });
        }

        // invalidate signatures of affected states
        for (auto currentState : statesWithInvalidBucketKey) {
            stateToBucketKey[currentState].reset();
        }
        statesWithInvalidBucketKey.clear();

        if (storm::utility::resources::isTerminate()) {
            // std::cout << "Performed " << iterations << " iterations of partition refinement before abort.\n";
            STORM_LOG_THROW(false, storm::exceptions::AbortException, "Aborted in bisimulation computation.");
            break;
        }
    }
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performEpsilonSignatureRefinementUsingCompleteLinkage(double epsilon) {
    std::cout << "Computing epsilon-bisimulation with epsilon: " << epsilon << std::endl;

    if constexpr (!storm::IsIntervalType<ValueType>) {
        STORM_LOG_THROW(true, storm::exceptions::IllegalFunctionCallException,
                        "Cannot compute epsilon-bisimulation on non-interval model!");
    }

    std::deque<typename Partition::Block> blocksQueue;
    storm::storage::Partition::BlockSet enqueuedBlocks;

    // enqueue all non-absorbing blocks
    this->partition.forEachBlock([&](auto const& block) {
        if (!this->absorbingBlocks.contains(block.front())) {
            blocksQueue.push_back(block);
            enqueuedBlocks.insert(block);
        }
    });

    const std::size_t numberOfStates = this->model.getTransitionMatrix().getRowCount();
    std::vector<ValueType> cachedTotalInterval(numberOfStates, storm::utility::zero<ValueType>());
    std::vector<std::unordered_map<storm::storage::sparse::state_type, ValueType>> cachedIntervalsToBlock(numberOfStates);

    // refinement loop
    uint_fast64_t iterations = 0;
    while (!blocksQueue.empty()) {
        ++iterations;

        auto block = blocksQueue.back();
        blocksQueue.pop_back();
        enqueuedBlocks.erase(block);

        refineBlockBasedOnEpsilonSignature(block, blocksQueue, enqueuedBlocks, epsilon);
    }

    std::cout << "Finished epsilon-refinement after " << iterations << " iterations." << std::endl;
    std::cout << "Size of final partition " << partition.getNumberOfBlocks() << "." << std::endl;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::performSignatureRefinement() {
    // Insert all blocks into the queue for refinement
    std::deque<typename Partition::Block> blocksQueue;

    storm::storage::Partition::BlockSet enqueuedSplitterBlocks;
    // Initially, add all current blocks to the queue
    this->partition.forEachBlock([&](auto const& block) {
        // TODO: maybe one has to handle the absorbing blocks differently for weak bisimulation, i.e., that they are still enqueued here and handled differently
        // while splitting based on the computed signatures
        if (!this->absorbingBlocks.contains(block.front())) {
            blocksQueue.push_back(block);
            enqueuedSplitterBlocks.insert(block);
        }
    });

    std::vector<size_t> stateToSignature(backwardTransitions.getColumnCount(), 0);
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
        auto wasSplit = this->partition.splitBlockByOrder(blockToRefine, [&stateToSignature]
                                                          (auto const& a, auto const& b) { return stateToSignature.at(a) < stateToSignature.at(b); });

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
                    for (auto &transition: backwardTransitions.getRow(state)) {
                        auto predecessorState = transition.getColumn();
                        auto predecessorBlock = partition.getBlockOfElement(predecessorState);

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
            // std::cout << "Performed " << iterations << " iterations of partition refinement before abort.\n";
            STORM_LOG_THROW(false, storm::exceptions::AbortException, "Aborted in bisimulation computation.");
            break;
        }
    }

    std::cout << "Finished refinement after " << iterations << " iterations." << std::endl;
    std::cout << "Attempt to split block failed " << noSplitCounter << " times." << std::endl;
}

template<typename ModelType>
std::size_t BisimulationDecomposition<ModelType>::computeStateSignatureHash(
    storm::storage::sparse::state_type state) const {
    // TODO: Interestingly, the boost hash function seems to create collisions, hence the quotient is not minimal per default.
    // TODO: When using --exact, the quotient gets calculated correctly. Thus, Boost does not guarantee a sufficient
    // TODO: precision for doubles
    return computeStateSignature(state, partition).computeHash();
    // TODO: Using --exact with the string based hash representation does not lead to a correct quotient -> investigate
    // return std::hash<std::string>{}(computeStateSignature(state, partition).toString());
}

template<typename ModelType>
storm::storage::bisimulation::Signature<typename ModelType::ValueType> BisimulationDecomposition<ModelType>::computeStateSignature(
    storm::storage::sparse::state_type state,
    storm::storage::bisimulation::Partition const& currentPartition) const {
    storm::storage::bisimulation::Signature<typename ModelType::ValueType> signature;

    for (auto entry : model.getTransitionMatrix().getRow(state)) {
        // std::cout << "Prob for state " << state << " to reach target state " << entry.getColumn() << ": " << entry.getValue() << std::endl;
        auto targetBlock = partition.getBlockOfElement(entry.getColumn()); // column marks the id of the target state
        signature.addBlockProbability(targetBlock.front(), entry.getValue());
    }

    return signature;
}

template<typename ModelType>
std::shared_ptr<ModelType> BisimulationDecomposition<ModelType>::getQuotient() const {
    STORM_LOG_THROW(this->quotient != nullptr, storm::exceptions::IllegalFunctionCallException,
                    "Unable to retrieve quotient model from bisimulation decomposition, because it was not built.");
    return this->quotient;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnRewards() {
    auto const& rewardModel = model.getUniqueRewardModel();
    if (rewardModel.hasStateRewards()) {
        this->splitInitialPartitionBasedOnRewards(rewardModel.getStateRewardVector());
    }
    if (rewardModel.hasStateActionRewards()) {
        if (model.isNondeterministicModel()) {
            std::vector<std::set<ValueType>> actionRewards;
            actionRewards.reserve(model.getNumberOfStates());
            for (storm::storage::sparse::state_type state = 0; state < model.getNumberOfStates(); ++state) {
                std::set<ValueType> rewardsAtState;
                for (auto choice = model.getTransitionMatrix().getRowGroupIndices()[state];
                     choice < model.getTransitionMatrix().getRowGroupIndices()[state + 1]; ++choice) {
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
    partition.forEachBlock([this, &rewardVector](auto const& block) {
        partition.splitBlockByOrder(block, [&rewardVector](auto const& a, auto const& b) { return rewardVector[a] < rewardVector[b]; });
    });
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::splitInitialPartitionBasedOnActionRewards(std::vector<std::set<ValueType>> const& actionRewards) {
    partition.forEachBlock([this, &actionRewards](auto const& block) {
        partition.splitBlockByOrder(block, [&actionRewards](auto const& a, auto const& b) { return actionRewards[a] < actionRewards[b]; });
    });
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::initializeLabelBasedPartition() {
    for (auto const& label : options.respectedAtomicPropositions.get()) {
        if (label == "init") {
            continue;
        }

        auto labelledStates = model.getStates(label);
        partition.forEachBlock([this, &labelledStates](auto const& block) {
            partition.splitBlockByPredicate(block, [&labelledStates](auto const& e) { return labelledStates.get(e); });
        });
    }

    // If the model has state rewards, we need to consider them, because otherwise reward properties are not
    // preserved.
    if (options.getKeepRewards() && model.hasRewardModel()) {
        // TODO: Check if this is implemented correctly
        this->splitInitialPartitionBasedOnRewards();
    }

    std::cout << "Number of blocks after initial partitioning by labels: " << partition.getNumberOfBlocks() << std::endl;
}

template<typename ModelType>
void BisimulationDecomposition<ModelType>::initializeMeasureDrivenPartition() {
    std::pair<storm::storage::BitVector, storm::storage::BitVector> statesWithProbability01 = this->getStatesWithProbability01();

    boost::optional<storm::storage::sparse::state_type> representativePsiState;
    if (!options.psiStates.get().empty()) {
        representativePsiState = *options.psiStates.get().begin();
    }

    auto prob0States = statesWithProbability01.first;
    auto prob1States = statesWithProbability01.second;

    // TODO: These blocks should be marked as absorbing
    auto initialBlock = partition.getBlockOfElement(0);
    if (!prob0States.empty()) {
        auto [notProb0Block, prob0Block] = partition.splitBlockByPredicate(initialBlock, [&prob0States]
                                                                           (auto const& state) { return prob0States.get(state); });
        absorbingBlocks.emplace(prob0Block.front(), prob0Block.front());

        if (!prob1States.empty()) {
            auto [notProb1Block, prob1Block] = partition.splitBlockByPredicate(notProb0Block, [&prob1States]
                                                                               (auto const& state) { return prob1States.get(state); });
            absorbingBlocks.emplace(prob1Block.front(), representativePsiState.get());
        }
    } else if (!prob1States.empty()) {
        auto [notProb1Block, prob1Block] = partition.splitBlockByPredicate(initialBlock, [&prob1States]
                                                                           (auto const& state) { return prob1States.get(state); });
        absorbingBlocks.emplace(prob1Block.front(), representativePsiState.get());

        if (!prob0States.empty()) {
            auto [notProb0Block, prob0Block] = partition.splitBlockByPredicate(notProb1Block, [&prob0States]
                                                                               (auto const& state) { return prob0States.get(state); });
            absorbingBlocks.emplace(prob0Block.front(), prob0Block.front());
        }
    }

    // If the model has state rewards, we need to consider them, because otherwise reward properties are not
    // preserved.
    if (options.getKeepRewards() && model.hasRewardModel()) {
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

#ifdef STORM_HAVE_CARL
template class BisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalNumber>>;
template class BisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalNumber>>;
template class BisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalNumber>>;

template class BisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalFunction>>;
template class BisimulationDecomposition<storm::models::sparse::Ctmc<storm::RationalFunction>>;
template class BisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalFunction>>;

template class storm::storage::BisimulationDecomposition<
    storm::models::sparse::Dtmc<carl::Interval<double>>>;
template class storm::storage::BisimulationDecomposition<
    storm::models::sparse::Ctmc<carl::Interval<double>>>;
template class storm::storage::BisimulationDecomposition<
    storm::models::sparse::Mdp<carl::Interval<double>>>;

#endif
}  // namespace storage
}  // namespace storm

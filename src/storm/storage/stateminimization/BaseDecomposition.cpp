#include "storm/storage/stateminimization/BaseDecomposition.h"

namespace storm {
namespace storage {

template<typename ModelType>
BaseDecomposition<ModelType>::BaseOptions::BaseOptions()
    : measureDrivenInitialPartition(false),
      phiStates(),
      psiStates(),
      respectedAtomicPropositions(),
      buildQuotient(true),
      keepRewards(false),
      bounded(false),
      discounted(false) {
    // Intentionally left empty.
}

template<typename ModelType>
BaseDecomposition<ModelType>::BaseOptions::BaseOptions(const ModelType &model, const storm::logic::Formula &formula) : BaseOptions() {
    this->preserveSingleFormula(model, formula);
}

template<typename ModelType>
BaseDecomposition<ModelType>::BaseOptions::BaseOptions(const ModelType &model, const std::vector<std::shared_ptr<const storm::logic::Formula>> &formulas)
    : BaseOptions() {
    if (formulas.empty()) {
        this->respectedAtomicPropositions = model.getStateLabeling().getLabels();
        this->keepRewards = true;
    }
    if (formulas.size() == 1) {
        this->preserveSingleFormula(model, *formulas.front());
    } else {
        for (auto const &formula : formulas) {
            preserveFormula(*formula);
        }
    }
}

template<typename ModelType>
void BaseDecomposition<ModelType>::BaseOptions::preserveSingleFormula(ModelType const &model, storm::logic::Formula const &formula) {
    // Retrieve information about formula.
    storm::logic::FormulaInformation info = formula.info();

    keepRewards = info.containsRewardOperator() || info.containsRewardBoundedFormula();

    // We need to preserve bounded properties iff the formula contains a bounded until or a next subformula.
    bounded = info.containsBoundedUntilFormula() || info.containsNextFormula() || info.containsCumulativeRewardFormula();

    // We need to preserve discounting iff the formula contains a discounted subformula.
    discounted = info.containsDiscountFormula();

    // Compute the relevant labels and expressions.
    this->addToRespectedAtomicPropositions(formula.getAtomicExpressionFormulas(), formula.getAtomicLabelFormulas());

    // Check whether measure driven initial partition is possible and, if so, set it.
    this->checkAndSetMeasureDrivenInitialPartition(model, formula);
}

// TODO: Is it clean to have this in the base class?
// TODO: This might not work for interval models, see https://github.com/stormchecker/storm/issues/899.
template<typename ModelType>
void BaseDecomposition<ModelType>::BaseOptions::checkAndSetMeasureDrivenInitialPartition(ModelType const &model, storm::logic::Formula const &formula) {
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

        using SolutionType = storm::IntervalBaseType<ValueType>;
        phiStates = phiStatesCheckResult->template asExplicitQualitativeCheckResult<SolutionType>().getTruthValuesVector();
        psiStates = psiStatesCheckResult->template asExplicitQualitativeCheckResult<SolutionType>().getTruthValuesVector();
    } else {
        optimalityType.reset();
    }
}

template<typename ModelType>
void BaseDecomposition<ModelType>::BaseOptions::addToRespectedAtomicPropositions(
    std::vector<std::shared_ptr<storm::logic::AtomicExpressionFormula const>> const &expressions,
    std::vector<std::shared_ptr<storm::logic::AtomicLabelFormula const>> const &labels) {
    std::set<std::string> labelsToRespect;
    for (auto const &labelFormula : labels) {
        labelsToRespect.insert(labelFormula->getLabel());
    }
    for (auto const &expressionFormula : expressions) {
        labelsToRespect.insert(expressionFormula->toString());
    }
    if (!respectedAtomicPropositions) {
        respectedAtomicPropositions = labelsToRespect;
    } else {
        respectedAtomicPropositions.value().insert(labelsToRespect.begin(), labelsToRespect.end());
    }
}

template<typename ModelType>
void BaseDecomposition<ModelType>::BaseOptions::preserveFormula(storm::logic::Formula const &formula) {
    // Disable the measure driven initial partition.
    measureDrivenInitialPartition = false;
    phiStates.reset();
    psiStates.reset();

    // Retrieve information about formula.
    storm::logic::FormulaInformation info = formula.info();

    // Preserve rewards if necessary.
    keepRewards = keepRewards || info.containsRewardOperator() || info.containsRewardBoundedFormula();

    // Preserve bounded properties if necessary.
    bounded = bounded || (info.containsBoundedUntilFormula() || info.containsNextFormula() || info.containsCumulativeRewardFormula());

    // Preserve discounted properties if necessary.
    discounted = discounted || info.containsDiscountFormula();

    // Compute the relevant labels and expressions.
    this->addToRespectedAtomicPropositions(formula.getAtomicExpressionFormulas(), formula.getAtomicLabelFormulas());
}

template<typename ModelType>
BaseDecomposition<ModelType>::BaseDecomposition(const ModelType &model, const storm::storage::SparseMatrix<ValueType> &backwardTransitions)
    : model(model), backwardTransitions(backwardTransitions), partition(model.getNumberOfStates()), quotient(nullptr), absorbingBlocks() {}

template<typename ModelType>
void BaseDecomposition<ModelType>::computeDecomposition() {
    auto totalStart = Clock::now();

    auto initialPartitionTime = this->measure([&]() { this->computeInitialPartition(); });

    STORM_LOG_WARN_COND(this->partition.getNumberOfBlocks() > 1, "Initial partition consists only of a single block.");

    auto refinementTime = this->measure([&]() { this->refinePartition(); });

    auto quotientBuildTime = this->measure([&]() {
        if (this->shouldBuildQuotient()) {
            this->buildQuotientFromPartition();
        }
    });

    auto totalTime = Clock::now() - totalStart;

    if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isShowStatisticsSet()) {
        this->printStatistics(initialPartitionTime, refinementTime, quotientBuildTime, totalTime);
    }
}

template<typename ModelType>
std::shared_ptr<ModelType> BaseDecomposition<ModelType>::getQuotient() const {
    STORM_LOG_THROW(this->quotient != nullptr, storm::exceptions::IllegalFunctionCallException,
                    "Unable to retrieve quotient model from decomposition, because it was not built.");
    return this->quotient;
}

template<typename ModelType>
void BaseDecomposition<ModelType>::printStatistics(Duration initialPartitionTime, Duration refinementTime, Duration quotientBuildTime,
                                                   Duration totalTime) const {
    auto toMilliseconds = [](Duration duration) { return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count(); };

    std::cout << "\nTime breakdown:\n";
    std::cout << "    * time for initial partition: " << toMilliseconds(initialPartitionTime) << "ms\n";
    std::cout << "    * time for partitioning: " << toMilliseconds(refinementTime) << "ms\n";
    std::cout << "    * time for building quotient: " << toMilliseconds(quotientBuildTime) << "ms\n";
    std::cout << "------------------------------------------\n";
    std::cout << "    * total time: " << toMilliseconds(totalTime) << "ms\n\n";
}

template class BaseDecomposition<storm::models::sparse::Dtmc<double>>;
template class BaseDecomposition<storm::models::sparse::Ctmc<double>>;
template class BaseDecomposition<storm::models::sparse::Mdp<double>>;

template class BaseDecomposition<storm::models::sparse::Dtmc<storm::RationalNumber>>;
template class BaseDecomposition<storm::models::sparse::Ctmc<storm::RationalNumber>>;
template class BaseDecomposition<storm::models::sparse::Mdp<storm::RationalNumber>>;

template class BaseDecomposition<storm::models::sparse::Dtmc<storm::RationalFunction>>;
template class BaseDecomposition<storm::models::sparse::Ctmc<storm::RationalFunction>>;
template class BaseDecomposition<storm::models::sparse::Mdp<storm::RationalFunction>>;

template class BaseDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;
template class BaseDecomposition<storm::models::sparse::Ctmc<storm::Interval>>;
template class BaseDecomposition<storm::models::sparse::Mdp<storm::Interval>>;

template class BaseDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>;
template class BaseDecomposition<storm::models::sparse::Ctmc<storm::RationalInterval>>;
template class BaseDecomposition<storm::models::sparse::Mdp<storm::RationalInterval>>;

}  // namespace storage
}  // namespace storm
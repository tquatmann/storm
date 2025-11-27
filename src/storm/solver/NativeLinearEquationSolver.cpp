#include "storm/solver/NativeLinearEquationSolver.h"

#include <limits>

#include "storm/environment/solver/NativeSolverEnvironment.h"
#include "storm/environment/solver/OviSolverEnvironment.h"

#include "storm/exceptions/InvalidEnvironmentException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnmetRequirementException.h"
#include "storm/solver/helper/GuessingValueIterationHelper.h"
#include "storm/solver/helper/IntervalterationHelper.h"
#include "storm/solver/helper/OptimisticValueIterationHelper.h"
#include "storm/solver/helper/RationalSearchHelper.h"
#include "storm/solver/helper/SoundValueIterationHelper.h"
#include "storm/solver/helper/ValueIterationHelper.h"
#include "storm/solver/multiplier/Multiplier.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/NumberTraits.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/constants.h"
#include "storm/utility/vector.h"

namespace storm {
namespace solver {

template<typename ValueType>
NativeLinearEquationSolver<ValueType>::NativeLinearEquationSolver() : localA(nullptr), A(nullptr) {
    // Intentionally left empty.
}

template<typename ValueType>
NativeLinearEquationSolver<ValueType>::NativeLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) : localA(nullptr), A(nullptr) {
    this->setMatrix(A);
}

template<typename ValueType>
NativeLinearEquationSolver<ValueType>::NativeLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) : localA(nullptr), A(nullptr) {
    this->setMatrix(std::move(A));
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
    localA.reset();
    this->A = &A;
    clearCache();
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
    localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
    this->A = localA.get();
    clearCache();
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::setUpViOperator() const {
    if (!viOperator) {
        viOperator = std::make_shared<helper::ValueIterationOperator<ValueType, true>>();
        viOperator->setMatrixBackwards(*this->A);
    }
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsSOR(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b,
                                                              ValueType const& omega) const {
    STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with NativeLinearEquationSolver (Gauss-Seidel, SOR omega = " << omega << ")");

    if (!this->cachedRowVector) {
        this->cachedRowVector = std::make_unique<std::vector<ValueType>>(getMatrixRowCount());
    }

    ValueType precision = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    uint64_t maxIter = env.solver().native().getMaximalNumberOfIterations();
    bool relative = env.solver().native().getRelativeTerminationCriterion();

    // Set up additional environment variables.
    uint_fast64_t iterations = 0;
    SolverStatus status = SolverStatus::InProgress;

    this->startMeasureProgress();
    while (status == SolverStatus::InProgress && iterations < maxIter) {
        A->performSuccessiveOverRelaxationStep(omega, x, b);

        // Now check if the process already converged within our precision.
        if (storm::utility::vector::equalModuloPrecision<ValueType>(*this->cachedRowVector, x, precision, relative)) {
            status = SolverStatus::Converged;
        }
        // If we did not yet converge, we need to backup the contents of x.
        if (status != SolverStatus::Converged) {
            *this->cachedRowVector = x;
        }

        // Potentially show progress.
        this->showProgressIterative(iterations);

        // Increase iteration count so we can abort if convergence is too slow.
        ++iterations;

        status = this->updateStatus(status, x, SolverGuarantee::None, iterations, maxIter);
    }

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    this->reportStatus(status, iterations);

    return status == SolverStatus::Converged;
}

template<typename ValueType>
NativeLinearEquationSolver<ValueType>::JacobiDecomposition::JacobiDecomposition(Environment const& env, storm::storage::SparseMatrix<ValueType> const& A) {
    auto decomposition = A.getJacobiDecomposition();
    this->LUMatrix = std::move(decomposition.first);
    this->DVector = std::move(decomposition.second);
    this->multiplier = storm::solver::MultiplierFactory<ValueType>().create(env, this->LUMatrix);
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsJacobi(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
    STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with NativeLinearEquationSolver (Jacobi)");

    if (!this->cachedRowVector) {
        this->cachedRowVector = std::make_unique<std::vector<ValueType>>(getMatrixRowCount());
    }

    // Get a Jacobi decomposition of the matrix A.
    if (!jacobiDecomposition) {
        jacobiDecomposition = std::make_unique<JacobiDecomposition>(env, *A);
    }

    ValueType precision = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    uint64_t maxIter = env.solver().native().getMaximalNumberOfIterations();
    bool relative = env.solver().native().getRelativeTerminationCriterion();

    std::vector<ValueType>* currentX = &x;
    std::vector<ValueType>* nextX = this->cachedRowVector.get();

    // Set up additional environment variables.
    uint_fast64_t iterations = 0;
    SolverStatus status = SolverStatus::InProgress;

    this->startMeasureProgress();
    while (status == SolverStatus::InProgress && iterations < maxIter) {
        // Compute D^-1 * (b - LU * x) and store result in nextX.
        jacobiDecomposition->multiplier->multiply(env, *currentX, nullptr, *nextX);
        storm::utility::vector::subtractVectors(b, *nextX, *nextX);
        storm::utility::vector::multiplyVectorsPointwise(jacobiDecomposition->DVector, *nextX, *nextX);

        // Now check if the process already converged within our precision.
        if (storm::utility::vector::equalModuloPrecision<ValueType>(*currentX, *nextX, precision, relative)) {
            status = SolverStatus::Converged;
        }
        // Swap the two pointers as a preparation for the next iteration.
        std::swap(nextX, currentX);

        // Potentially show progress.
        this->showProgressIterative(iterations);

        // Increase iteration count so we can abort if convergence is too slow.
        ++iterations;

        status = this->updateStatus(status, *currentX, SolverGuarantee::None, iterations, maxIter);
    }

    // If the last iteration did not write to the original x we have to swap the contents, because the
    // output has to be written to the input parameter x.
    if (currentX == this->cachedRowVector.get()) {
        std::swap(x, *currentX);
    }

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    this->reportStatus(status, iterations);

    return status == SolverStatus::Converged;
}

template<typename ValueType>
NativeLinearEquationSolver<ValueType>::WalkerChaeData::WalkerChaeData(Environment const& env, storm::storage::SparseMatrix<ValueType> const& originalMatrix,
                                                                      std::vector<ValueType> const& originalB)
    : t(storm::utility::convertNumber<ValueType>(1000.0)) {
    computeWalkerChaeMatrix(originalMatrix);
    computeNewB(originalB);
    precomputeAuxiliaryData();
    multiplier = storm::solver::MultiplierFactory<ValueType>().create(env, this->matrix);
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::WalkerChaeData::computeWalkerChaeMatrix(storm::storage::SparseMatrix<ValueType> const& originalMatrix) {
    storm::storage::BitVector columnsWithNegativeEntries(originalMatrix.getColumnCount());
    ValueType zero = storm::utility::zero<ValueType>();
    for (auto const& e : originalMatrix) {
        if (e.getValue() < zero) {
            columnsWithNegativeEntries.set(e.getColumn());
        }
    }
    std::vector<uint64_t> columnsWithNegativeEntriesBefore = columnsWithNegativeEntries.getNumberOfSetBitsBeforeIndices();

    // We now build an extended equation system matrix that only has non-negative coefficients.
    storm::storage::SparseMatrixBuilder<ValueType> builder;

    uint64_t row = 0;
    for (; row < originalMatrix.getRowCount(); ++row) {
        for (auto const& entry : originalMatrix.getRow(row)) {
            if (entry.getValue() < zero) {
                builder.addNextValue(row, originalMatrix.getRowCount() + columnsWithNegativeEntriesBefore[entry.getColumn()], -entry.getValue());
            } else {
                builder.addNextValue(row, entry.getColumn(), entry.getValue());
            }
        }
    }
    ValueType one = storm::utility::one<ValueType>();
    for (auto column : columnsWithNegativeEntries) {
        builder.addNextValue(row, column, one);
        builder.addNextValue(row, originalMatrix.getRowCount() + columnsWithNegativeEntriesBefore[column], one);
        ++row;
    }

    matrix = builder.build();
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::WalkerChaeData::computeNewB(std::vector<ValueType> const& originalB) {
    b = std::vector<ValueType>(originalB);
    b.resize(matrix.getRowCount());
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::WalkerChaeData::precomputeAuxiliaryData() {
    columnSums = std::vector<ValueType>(matrix.getColumnCount());
    for (auto const& e : matrix) {
        columnSums[e.getColumn()] += e.getValue();
    }

    newX.resize(matrix.getRowCount());
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsWalkerChae(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
    STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with NativeLinearEquationSolver (WalkerChae)");

    // (1) Compute an equivalent equation system that has only non-negative coefficients.
    if (!walkerChaeData) {
        walkerChaeData = std::make_unique<WalkerChaeData>(env, *this->A, b);
    }

    // (2) Enlarge the vectors x and b to account for additional variables.
    x.resize(walkerChaeData->matrix.getRowCount());

    // Square the error bound, so we can use it to check for convergence. We take the squared error, because we
    // do not want to compute the root in the 2-norm computation.
    ValueType squaredErrorBound = storm::utility::pow(storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision()), 2);

    uint64_t maxIter = env.solver().native().getMaximalNumberOfIterations();

    // Set up references to the x-vectors used in the iteration loop.
    std::vector<ValueType>* currentX = &x;
    std::vector<ValueType>* nextX = &walkerChaeData->newX;

    std::vector<ValueType> tmp = walkerChaeData->matrix.getRowSumVector();
    storm::utility::vector::applyPointwise(tmp, walkerChaeData->b, walkerChaeData->b,
                                           [this](ValueType const& first, ValueType const& second) -> ValueType { return walkerChaeData->t * first + second; });

    // Add t to all entries of x.
    storm::utility::vector::applyPointwise(x, x, [this](ValueType const& value) -> ValueType { return value + walkerChaeData->t; });

    // Create a vector that always holds Ax.
    std::vector<ValueType> currentAx(x.size());
    walkerChaeData->multiplier->multiply(env, *currentX, nullptr, currentAx);

    // (3) Perform iterations until convergence.
    SolverStatus status = SolverStatus::InProgress;
    uint64_t iterations = 0;
    this->startMeasureProgress();
    while (status == SolverStatus::InProgress && iterations < maxIter) {
        // Perform one Walker-Chae step.
        walkerChaeData->matrix.performWalkerChaeStep(*currentX, walkerChaeData->columnSums, walkerChaeData->b, currentAx, *nextX);

        // Compute new Ax.
        walkerChaeData->multiplier->multiply(env, *nextX, nullptr, currentAx);

        // Check for convergence.
        if (storm::utility::vector::computeSquaredNorm2Difference(currentAx, walkerChaeData->b) <= squaredErrorBound) {
            status = SolverStatus::Converged;
        }

        // Swap the x vectors for the next iteration.
        std::swap(currentX, nextX);

        // Potentially show progress.
        this->showProgressIterative(iterations);

        // Increase iteration count so we can abort if convergence is too slow.
        ++iterations;

        if (storm::utility::resources::isTerminate()) {
            status = SolverStatus::Aborted;
        }
    }

    // If the last iteration did not write to the original x we have to swap the contents, because the
    // output has to be written to the input parameter x.
    if (currentX == &walkerChaeData->newX) {
        std::swap(x, *currentX);
    }

    // Resize the solution to the right size.
    x.resize(this->A->getRowCount());

    // Finalize solution vector.
    storm::utility::vector::applyPointwise(x, x, [this](ValueType const& value) -> ValueType { return value - walkerChaeData->t; });

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    this->reportStatus(status, iterations);

    return status == SolverStatus::Converged;
}

template<typename ValueType>
typename NativeLinearEquationSolver<ValueType>::PowerIterationResult NativeLinearEquationSolver<ValueType>::performPowerIteration(
    Environment const& env, std::vector<ValueType>*& currentX, std::vector<ValueType>*& newX, std::vector<ValueType> const& b, ValueType const& precision,
    bool relative, SolverGuarantee const& guarantee, uint64_t currentIterations, uint64_t maxIterations,
    storm::solver::MultiplicationStyle const& multiplicationStyle) const {
    bool useGaussSeidelMultiplication = multiplicationStyle == storm::solver::MultiplicationStyle::GaussSeidel;

    uint64_t iterations = currentIterations;
    SolverStatus status = this->terminateNow(*currentX, guarantee) ? SolverStatus::TerminatedEarly : SolverStatus::InProgress;
    while (status == SolverStatus::InProgress && iterations < maxIterations) {
        if (useGaussSeidelMultiplication) {
            *newX = *currentX;
            this->multiplier->multiplyGaussSeidel(env, *newX, &b);
        } else {
            this->multiplier->multiply(env, *currentX, &b, *newX);
        }

        // Check for convergence.
        if (storm::utility::vector::equalModuloPrecision<ValueType>(*currentX, *newX, precision, relative)) {
            status = SolverStatus::Converged;
        }

        // Check for termination.
        std::swap(currentX, newX);
        ++iterations;

        status = this->updateStatus(status, *currentX, guarantee, iterations, maxIterations);

        // Potentially show progress.
        this->showProgressIterative(iterations);
    }

    return PowerIterationResult(iterations - currentIterations, status);
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsPower(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
    STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with NativeLinearEquationSolver (Power)");
    // Prepare the solution vectors.
    setUpViOperator();

    SolverGuarantee guarantee = SolverGuarantee::None;
    if (this->hasCustomTerminationCondition()) {
        if (this->getTerminationCondition().requiresGuarantee(SolverGuarantee::LessOrEqual) && this->hasLowerBound()) {
            this->createLowerBoundsVector(x);
            guarantee = SolverGuarantee::LessOrEqual;
        } else if (this->getTerminationCondition().requiresGuarantee(SolverGuarantee::GreaterOrEqual) && this->hasUpperBound()) {
            this->createUpperBoundsVector(x);
            guarantee = SolverGuarantee::GreaterOrEqual;
        }
    }

    storm::solver::helper::ValueIterationHelper<ValueType, true> viHelper(viOperator);
    uint64_t numIterations{0};
    auto viCallback = [&](SolverStatus const& current) {
        this->showProgressIterative(numIterations);
        return this->updateStatus(current, x, guarantee, numIterations, env.solver().native().getMaximalNumberOfIterations());
    };
    this->startMeasureProgress();
    auto status = viHelper.VI(x, b, numIterations, env.solver().native().getRelativeTerminationCriterion(),
                              storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision()), {}, viCallback,
                              env.solver().native().getPowerMethodMultiplicationStyle());

    this->reportStatus(status, numIterations);

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
void preserveOldRelevantValues(std::vector<ValueType> const& allValues, storm::storage::BitVector const& relevantValues, std::vector<ValueType>& oldValues) {
    storm::utility::vector::selectVectorValues(oldValues, relevantValues, allValues);
}

template<typename ValueType>
ValueType computeMaxAbsDiff(std::vector<ValueType> const& allValues, storm::storage::BitVector const& relevantValues, std::vector<ValueType> const& oldValues) {
    ValueType result = storm::utility::zero<ValueType>();
    auto oldValueIt = oldValues.begin();
    for (auto value : relevantValues) {
        result = storm::utility::max<ValueType>(result, storm::utility::abs<ValueType>(allValues[value] - *oldValueIt));
    }
    return result;
}

template<typename ValueType>
ValueType computeMaxAbsDiff(std::vector<ValueType> const& allOldValues, std::vector<ValueType> const& allNewValues,
                            storm::storage::BitVector const& relevantValues) {
    ValueType result = storm::utility::zero<ValueType>();
    for (auto value : relevantValues) {
        result = storm::utility::max<ValueType>(result, storm::utility::abs<ValueType>(allNewValues[value] - allOldValues[value]));
    }
    return result;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsIntervalIteration(Environment const& env, std::vector<ValueType>& xLower,
                                                                            std::vector<ValueType>& xUpper, std::vector<ValueType> const& bLower,
                                                                            std::vector<ValueType> const& bUpper) const {
    STORM_LOG_THROW(this->hasLowerBound(), storm::exceptions::UnmetRequirementException, "Solver requires lower bound, but none was given.");
    STORM_LOG_THROW(this->hasUpperBound(), storm::exceptions::UnmetRequirementException, "Solver requires upper bound, but none was given.");
    STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with NativeLinearEquationSolver (IntervalIteration)");
    setUpViOperator();
    helper::IntervalIterationHelper<ValueType, true> iiHelper(viOperator);
    auto prec = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    auto lowerBoundsCallback = [&](std::vector<ValueType>& vector) { this->createLowerBoundsVector(vector); };
    auto upperBoundsCallback = [&](std::vector<ValueType>& vector) { this->createUpperBoundsVector(vector); };

    uint64_t numIterations{0};
    auto iiCallback = [&](helper::IIData<ValueType> const& data) {
        this->showProgressIterative(numIterations);
        bool terminateEarly = this->hasCustomTerminationCondition() && this->getTerminationCondition().terminateNow(data.x, SolverGuarantee::LessOrEqual) &&
                              this->getTerminationCondition().terminateNow(data.y, SolverGuarantee::GreaterOrEqual);
        return this->updateStatus(data.status, terminateEarly, numIterations, env.solver().native().getMaximalNumberOfIterations());
    };
    std::optional<storm::storage::BitVector> optionalRelevantValues;
    if (this->hasRelevantValues()) {
        optionalRelevantValues = this->getRelevantValues();
    }
    this->startMeasureProgress();
    SolverStatus status = SolverStatus::InProgress;
    if (&xLower != &xUpper) {
        std::pair<std::vector<ValueType>, std::vector<ValueType>> xy{std::move(xLower), std::move(xUpper)};
        lowerBoundsCallback(xy.first);
        upperBoundsCallback(xy.second);
        if (&bLower != &bUpper) {
            std::pair<std::vector<ValueType> const&, std::vector<ValueType> const&> b{bLower, bUpper};
            status = iiHelper.II(xy, b, numIterations, env.solver().native().getRelativeTerminationCriterion(), prec, {}, iiCallback, optionalRelevantValues);
        } else {
            status =
                iiHelper.II(xy, bLower, numIterations, env.solver().native().getRelativeTerminationCriterion(), prec, {}, iiCallback, optionalRelevantValues);
        }
        xLower = std::move(xy.first);
        xUpper = std::move(xy.second);
    } else {
        status = iiHelper.II(xLower, bLower, numIterations, env.solver().native().getRelativeTerminationCriterion(), prec, lowerBoundsCallback,
                             upperBoundsCallback, {}, iiCallback, optionalRelevantValues);
    }
    this->reportStatus(status, numIterations);

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsSoundValueIteration(Environment const& env, std::vector<ValueType>& xLower,
                                                                              std::vector<ValueType>& xUpper, std::vector<ValueType> const& bLower,
                                                                              std::vector<ValueType> const& bUpper) const {
    STORM_LOG_THROW(&bLower == &bUpper, storm::exceptions::NotSupportedException,
                    "Optimistic value iteration with different lower/upper b-vectors is not supported.");
    auto const& b = bLower;

    // Prepare the solution vectors and the helper.
    STORM_LOG_ASSERT(xLower.size() == this->A->getRowGroupCount(), "Unexpected size of xLower");
    STORM_LOG_ASSERT(xUpper.size() == this->A->getRowGroupCount(), "Unexpected size of xUpper");

    std::optional<ValueType> lowerBound, upperBound;
    if (this->hasLowerBound()) {
        lowerBound = this->getLowerBound(true);
    }
    if (this->hasUpperBound()) {
        upperBound = this->getUpperBound(true);
    }

    setUpViOperator();

    auto precision = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    uint64_t numIterations{0};
    auto sviCallback = [&](typename helper::SoundValueIterationHelper<ValueType, true>::SVIData const& current) {
        this->showProgressIterative(numIterations);
        return this->updateStatus(current.status,
                                  this->hasCustomTerminationCondition() && current.checkCustomTerminationCondition(this->getTerminationCondition()),
                                  numIterations, env.solver().native().getMaximalNumberOfIterations());
    };
    std::optional<storm::storage::BitVector> optionalRelevantValues;
    if (this->hasRelevantValues()) {
        optionalRelevantValues = this->getRelevantValues();
    }
    this->startMeasureProgress();
    helper::SoundValueIterationHelper<ValueType, true> sviHelper(viOperator);
    SolverStatus status = SolverStatus::InProgress;
    if (&xLower != &xUpper) {
        std::pair<std::vector<ValueType>, std::vector<ValueType>> xy{std::move(xLower), std::move(xUpper)};
        auto sviData = sviHelper.SVI(xy, b, numIterations, env.solver().native().getRelativeTerminationCriterion(), precision, {}, lowerBound, upperBound,
                                     sviCallback, optionalRelevantValues);
        sviData.trySetLowerUpper(xy.first, xy.second);
        xLower = std::move(xy.first);
        xUpper = std::move(xy.second);
        status = sviData.status;
    } else {
        status = sviHelper.SVI(xLower, b, numIterations, env.solver().native().getRelativeTerminationCriterion(), precision, {}, lowerBound, upperBound,
                               sviCallback, optionalRelevantValues);
    }

    this->reportStatus(status, numIterations);

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsOptimisticValueIteration(Environment const& env, std::vector<ValueType>& xLower,
                                                                                   std::vector<ValueType>& xUpper, std::vector<ValueType> const& bLower,
                                                                                   std::vector<ValueType> const& bUpper) const {
    STORM_LOG_THROW(&bLower == &bUpper, storm::exceptions::NotSupportedException,
                    "Optimistic value iteration with different lower/upper b-vectors is not supported.");
    auto const& b = bLower;

    if (!storm::utility::vector::hasNonZeroEntry(b)) {
        // If all entries are zero, OVI might run in an endless loop. However, the result is easy in this case.
        xLower.assign(xLower.size(), storm::utility::zero<ValueType>());
        if (&xLower != &xUpper) {
            xUpper = xLower;
        }
        return true;
    }

    setUpViOperator();

    helper::OptimisticValueIterationHelper<ValueType, true> oviHelper(viOperator);
    auto const prec = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    auto const guessingFactor =
        storm::utility::convertNumber<ValueType>(env.solver().ovi().getUpperBoundGuessingFactor().value_or(env.solver().native().getPrecision()));
    std::optional<ValueType> lowerBound, upperBound;
    if (this->hasLowerBound()) {
        lowerBound = this->getLowerBound(true);
    }
    if (this->hasUpperBound()) {
        upperBound = this->getUpperBound(true);
    }
    uint64_t numIterations{0};
    auto oviCallback = [&](SolverStatus const& current, std::vector<ValueType> const& v) {
        this->showProgressIterative(numIterations);
        return this->updateStatus(current, v, SolverGuarantee::LessOrEqual, numIterations, env.solver().native().getMaximalNumberOfIterations());
    };
    this->createLowerBoundsVector(xLower);
    this->startMeasureProgress();
    SolverStatus status = SolverStatus::InProgress;
    if (&xLower != &xUpper) {
        std::pair<std::vector<ValueType>, std::vector<ValueType>> vu{std::move(xLower), std::move(xUpper)};
        status = oviHelper.OVI(vu, b, numIterations, env.solver().native().getRelativeTerminationCriterion(), prec, {}, guessingFactor, lowerBound, upperBound,
                               oviCallback);
        xLower = std::move(vu.first);
        xUpper = std::move(vu.second);
    } else {
        status = oviHelper.OVI(xLower, b, numIterations, env.solver().native().getRelativeTerminationCriterion(), prec, {}, guessingFactor, lowerBound,
                               upperBound, oviCallback);
    }
    this->reportStatus(status, numIterations);

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsGuessingValueIteration(const Environment& env, std::vector<ValueType>& xLower,
                                                                                 std::vector<ValueType>& xUpper, std::vector<ValueType> const& bLower,
                                                                                 std::vector<ValueType> const& bUpper) const {
    STORM_LOG_THROW(&bLower == &bUpper, storm::exceptions::NotSupportedException,
                    "Optimistic value iteration with different lower/upper b-vectors is not supported.");
    auto const& b = bLower;

    setUpViOperator();

    if (&xLower == &xUpper && !this->cachedRowVector) {
        this->cachedRowVector = std::make_unique<std::vector<ValueType>>(this->A->getRowCount());
    }
    auto& xUpperRef = (&xLower != &xUpper) ? xUpper : *this->cachedRowVector;

    storm::solver::helper::GuessingValueIterationHelper<ValueType, true> helper(viOperator, *this->A);

    uint64_t numIterations{0};
    auto gviCallback = [&](helper::GVIData<ValueType> const& data) {
        this->showProgressIterative(numIterations);
        bool terminateEarly = this->hasCustomTerminationCondition() && this->getTerminationCondition().terminateNow(data.x, SolverGuarantee::LessOrEqual) &&
                              this->getTerminationCondition().terminateNow(data.y, SolverGuarantee::GreaterOrEqual);
        return this->updateStatus(data.status, terminateEarly, numIterations, env.solver().native().getMaximalNumberOfIterations());
    };

    this->createLowerBoundsVector(xLower);
    this->createUpperBoundsVector(xUpperRef);

    ValueType precision = storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision());
    auto const two = storm::utility::convertNumber<ValueType>(2.0);
    if (&xLower == &xUpper) {
        // It suffices to have one vector that is close enough. This means that the gap between lower/upper bound can be twice as large
        precision *= two;
    }
    bool const relative = env.solver().native().getRelativeTerminationCriterion();
    STORM_LOG_THROW(!relative, storm::exceptions::NotSupportedException,
                    "Guessing value iteration currently does not support relative precision. Try to set an absolute precision criterion.");

    this->startMeasureProgress();
    auto status = helper.solveEquations(xLower, xUpperRef, b, numIterations, precision, {}, gviCallback);
    this->reportStatus(status, numIterations);

    if (&xLower == &xUpper) {
        // Set the averaged values
        storm::utility::vector::applyPointwise<ValueType, ValueType, ValueType>(
            xLower, xUpper, xLower, [&two](ValueType const& a, ValueType const& b) -> ValueType { return (a + b) / two; });
    }

    if (!this->isCachingEnabled()) {
        clearCache();
    }
    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::solveEquationsRationalSearch(Environment const& env, std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                         std::vector<ValueType> const& bLower, std::vector<ValueType> const& bUpper) const {
    // Set up two value iteration operators. One for exact and one for imprecise computations
    setUpViOperator();
    std::shared_ptr<helper::ValueIterationOperator<storm::RationalNumber, true>> exactOp;
    std::shared_ptr<helper::ValueIterationOperator<double, true>> impreciseOp;

    if constexpr (std::is_same_v<ValueType, storm::RationalNumber>) {
        exactOp = viOperator;
        impreciseOp = std::make_shared<helper::ValueIterationOperator<double, true>>();
        impreciseOp->setMatrixBackwards(this->A->template toValueType<double>());
    } else {
        impreciseOp = viOperator;
        exactOp = std::make_shared<helper::ValueIterationOperator<storm::RationalNumber, true>>();
        exactOp->setMatrixBackwards(this->A->template toValueType<storm::RationalNumber>());
    }

    storm::solver::helper::RationalSearchHelper<ValueType, storm::RationalNumber, double, true> rsHelper(exactOp, impreciseOp);
    uint64_t numIterations{0};
    auto applyRationalSearch = [&](auto& x, auto const& b) {
        auto rsCallback = [&](SolverStatus const& current) {
            this->showProgressIterative(numIterations);
            return this->updateStatus(current, x, SolverGuarantee::None, numIterations, env.solver().native().getMaximalNumberOfIterations());
        };
        return rsHelper.RS(x, b, numIterations, storm::utility::convertNumber<ValueType>(env.solver().native().getPrecision()), {}, rsCallback);
    };

    this->startMeasureProgress();
    auto status = applyRationalSearch(xLower, bLower);
    if ((&xLower != &xUpper) && (status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly)) {
        if (&bLower == &bUpper) {
            xUpper = xLower;  // no need to recompute, just copy over values
        } else {
            status = applyRationalSearch(xUpper, bUpper);
        }
    }
    this->reportStatus(status, numIterations);

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return status == SolverStatus::Converged || status == SolverStatus::TerminatedEarly;
}

template<typename ValueType>
NativeLinearEquationSolverMethod NativeLinearEquationSolver<ValueType>::getMethod(Environment const& env, bool isExactMode) const {
    // Adjust the method if none was specified and we want exact or sound computations
    auto method = env.solver().native().getMethod();

    if (isExactMode && method != NativeLinearEquationSolverMethod::RationalSearch) {
        if (env.solver().native().isMethodSetFromDefault()) {
            method = NativeLinearEquationSolverMethod::RationalSearch;
            STORM_LOG_INFO(
                "Selecting '" + toString(method) +
                "' as the solution technique to guarantee exact results. If you want to override this, please explicitly specify a different method.");
        } else {
            STORM_LOG_WARN("The selected solution method does not guarantee exact results.");
        }
    } else if (env.solver().isForceSoundness() && method != NativeLinearEquationSolverMethod::SoundValueIteration &&
               method != NativeLinearEquationSolverMethod::OptimisticValueIteration && method != NativeLinearEquationSolverMethod::IntervalIteration &&
               method != NativeLinearEquationSolverMethod::RationalSearch) {
        if (env.solver().native().isMethodSetFromDefault()) {
            method = NativeLinearEquationSolverMethod::OptimisticValueIteration;
            STORM_LOG_INFO(
                "Selecting '" + toString(method) +
                "' as the solution technique to guarantee sound results. If you want to override this, please explicitly specify a different method.");
        } else {
            STORM_LOG_WARN("The selected solution method does not guarantee sound results.");
        }
    }
    return method;
}

template<typename ValueType>
bool NativeLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                   std::vector<ValueType> const& bLower, std::vector<ValueType> const& bUpper) const {
    STORM_LOG_ASSERT(&xLower != &xUpper || &bLower == &bUpper, "solving with different lower and upper b-values requires different lower and upper solutions.");

    switch (getMethod(env, storm::NumberTraits<ValueType>::IsExact || env.solver().isForceExact())) {
        case NativeLinearEquationSolverMethod::SOR:
            STORM_LOG_ASSERT(&xLower == &xUpper, "The selected solution method can not be used to compute sound bounds.");
            return this->solveEquationsSOR(env, xLower, bLower, storm::utility::convertNumber<ValueType>(env.solver().native().getSorOmega()));
        case NativeLinearEquationSolverMethod::GaussSeidel:
            STORM_LOG_ASSERT(&xLower == &xUpper, "The selected solution method can not be used to compute sound bounds.");
            return this->solveEquationsSOR(env, xLower, bLower, storm::utility::one<ValueType>());
        case NativeLinearEquationSolverMethod::Jacobi:
            STORM_LOG_ASSERT(&xLower == &xUpper, "The selected solution method can not be used to compute sound bounds.");
            return this->solveEquationsJacobi(env, xLower, bLower);
        case NativeLinearEquationSolverMethod::WalkerChae:
            STORM_LOG_ASSERT(&xLower == &xUpper, "The selected solution method can not be used to compute sound bounds.");
            return this->solveEquationsWalkerChae(env, xLower, bLower);
        case NativeLinearEquationSolverMethod::Power:
            STORM_LOG_ASSERT(&xLower == &xUpper, "The selected solution method can not be used to compute sound bounds.");
            return this->solveEquationsPower(env, xLower, bLower);
        case NativeLinearEquationSolverMethod::SoundValueIteration:
            return this->solveEquationsSoundValueIteration(env, xLower, xUpper, bLower, bUpper);
        case NativeLinearEquationSolverMethod::OptimisticValueIteration:
            return this->solveEquationsOptimisticValueIteration(env, xLower, xUpper, bLower, bUpper);
        case NativeLinearEquationSolverMethod::GuessingValueIteration:
            return this->solveEquationsGuessingValueIteration(env, xLower, xUpper, bLower, bUpper);
        case NativeLinearEquationSolverMethod::IntervalIteration:
            return this->solveEquationsIntervalIteration(env, xLower, xUpper, bLower, bUpper);
        case NativeLinearEquationSolverMethod::RationalSearch:
            return this->solveEquationsRationalSearch(env, xLower, xUpper, bLower, bUpper);
    }
    STORM_LOG_THROW(false, storm::exceptions::InvalidEnvironmentException, "Unknown solving technique.");
    return false;
}

template<typename ValueType>
LinearEquationSolverProblemFormat NativeLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const& env) const {
    auto method = getMethod(env, storm::NumberTraits<ValueType>::IsExact || env.solver().isForceExact());
    if (method == NativeLinearEquationSolverMethod::Power || method == NativeLinearEquationSolverMethod::SoundValueIteration ||
        method == NativeLinearEquationSolverMethod::OptimisticValueIteration || method == NativeLinearEquationSolverMethod::RationalSearch ||
        method == NativeLinearEquationSolverMethod::IntervalIteration || method == NativeLinearEquationSolverMethod::GuessingValueIteration) {
        return LinearEquationSolverProblemFormat::FixedPointSystem;
    } else {
        return LinearEquationSolverProblemFormat::EquationSystem;
    }
}

template<typename ValueType>
LinearEquationSolverRequirements NativeLinearEquationSolver<ValueType>::getRequirements(Environment const& env) const {
    LinearEquationSolverRequirements requirements;
    auto method = getMethod(env, storm::NumberTraits<ValueType>::IsExact || env.solver().isForceExact());
    if (method == NativeLinearEquationSolverMethod::IntervalIteration || method == NativeLinearEquationSolverMethod::GuessingValueIteration) {
        requirements.requireBounds();
    } else if (method == NativeLinearEquationSolverMethod::RationalSearch || method == NativeLinearEquationSolverMethod::OptimisticValueIteration) {
        requirements.requireLowerBounds();
    } else if (method == NativeLinearEquationSolverMethod::SoundValueIteration) {
        requirements.requireBounds(false);
    }
    return requirements;
}

template<typename ValueType>
void NativeLinearEquationSolver<ValueType>::clearCache() const {
    jacobiDecomposition.reset();
    walkerChaeData.reset();
    multiplier.reset();
    viOperator.reset();
    LinearEquationSolver<ValueType>::clearCache();
}

template<typename ValueType>
uint64_t NativeLinearEquationSolver<ValueType>::getMatrixRowCount() const {
    return this->A->getRowCount();
}

template<typename ValueType>
uint64_t NativeLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
    return this->A->getColumnCount();
}

template<typename ValueType>
std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> NativeLinearEquationSolverFactory<ValueType>::create(Environment const&) const {
    return std::make_unique<storm::solver::NativeLinearEquationSolver<ValueType>>();
}

template<typename ValueType>
std::unique_ptr<LinearEquationSolverFactory<ValueType>> NativeLinearEquationSolverFactory<ValueType>::clone() const {
    return std::make_unique<NativeLinearEquationSolverFactory<ValueType>>(*this);
}

// Explicitly instantiate the linear equation solver.
template class NativeLinearEquationSolver<double>;
template class NativeLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
template class NativeLinearEquationSolver<storm::RationalNumber>;
template class NativeLinearEquationSolverFactory<storm::RationalNumber>;

#endif
}  // namespace solver
}  // namespace storm

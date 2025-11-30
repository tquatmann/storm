#include "storm/solver/TopologicalLinearEquationSolver.h"

#include "storm/environment/solver/TopologicalSolverEnvironment.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/exceptions/InvalidEnvironmentException.h"
#include "storm/exceptions/InvalidStateException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/ProgressMeasurement.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/constants.h"
#include "storm/utility/vector.h"

namespace storm {
namespace solver {

template<typename ValueType>
TopologicalLinearEquationSolver<ValueType>::TopologicalLinearEquationSolver() : localA(nullptr), A(nullptr) {
    // Intentionally left empty.
}

template<typename ValueType>
TopologicalLinearEquationSolver<ValueType>::TopologicalLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) : localA(nullptr), A(nullptr) {
    this->setMatrix(A);
}

template<typename ValueType>
TopologicalLinearEquationSolver<ValueType>::TopologicalLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) : localA(nullptr), A(nullptr) {
    this->setMatrix(std::move(A));
}

template<typename ValueType>
void TopologicalLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
    localA.reset();
    this->A = &A;
    clearCache();
}

template<typename ValueType>
void TopologicalLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
    localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
    this->A = localA.get();
    clearCache();
}

inline storm::Environment getEnvironmentForUnderlyingSolver(storm::Environment const& env) {
    storm::Environment subEnv(env);
    subEnv.solver().setLinearEquationSolverType(env.solver().topological().getUnderlyingEquationSolverType(),
                                                env.solver().topological().isUnderlyingEquationSolverTypeSetFromDefault());
    return subEnv;
}

template<typename ValueType>
bool TopologicalLinearEquationSolver<ValueType>::supportsSolveEquationsSoundBounds(Environment const& env, bool) const {
    auto solverEnv = getEnvironmentForUnderlyingSolver(env);
    return GeneralLinearEquationSolverFactory<ValueType>().create(solverEnv)->supportsSolveEquationsSoundBounds(solverEnv, true);
}

template<typename ValueType>
bool TopologicalLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                        std::vector<ValueType> const& bLower, std::vector<ValueType> const& bUpper) const {
    STORM_LOG_ASSERT(&xLower != &xUpper || &bLower == &bUpper, "Solving with different lower and upper b-values requires different lower and upper solutions.");
    STORM_LOG_ASSERT(xLower.size() == this->A->getRowGroupCount(), "Provided x-vector (lower) has invalid size.");
    STORM_LOG_ASSERT(xUpper.size() == this->A->getRowGroupCount(), "Provided x-vector (upper) has invalid size.");
    STORM_LOG_ASSERT(bLower.size() == this->A->getRowCount(), "Provided b-vector (lower) has invalid size.");
    STORM_LOG_ASSERT(bUpper.size() == this->A->getRowCount(), "Provided b-vector (upper) has invalid size.");

    // For sound computations we need to increase the precision in each SCC
    storm::Environment sccSolverEnvironment = getEnvironmentForUnderlyingSolver(env);
    auto const [precision, relative] = env.solver().getPrecisionOfLinearEquationSolver(env.solver().topological().getUnderlyingEquationSolverType());
    bool needAdaptPrecision = env.solver().isForceSoundness() && precision.is_initialized();

    if (!this->sortedSccDecomposition || (needAdaptPrecision && !this->longestSccChainSize)) {
        STORM_LOG_TRACE("Creating SCC decomposition.");
        storm::utility::Stopwatch sccSw(true);
        createSortedSccDecomposition(needAdaptPrecision);
        sccSw.stop();
        STORM_LOG_INFO("SCC decomposition computed in "
                       << sccSw << ". Found " << this->sortedSccDecomposition->size() << " SCC(s) containing a total of " << xLower.size()
                       << " states. Average SCC size is "
                       << static_cast<double>(this->getMatrixRowCount()) / static_cast<double>(this->sortedSccDecomposition->size()) << ".");
    }

    // We do not need to adapt the precision if all SCCs are trivial (i.e., the system is acyclic)
    needAdaptPrecision = needAdaptPrecision && (this->sortedSccDecomposition->size() != this->getMatrixRowCount());
    auto setPrecisionDivisor = [&sccSolverEnvironment, &precision](uint64_t divisor) {
        sccSolverEnvironment.solver().setLinearEquationSolverPrecision(
            static_cast<storm::RationalNumber>((*precision) / storm::utility::convertNumber<storm::RationalNumber>(divisor)));
    };
    if (this->longestSccChainSize) {
        STORM_LOG_INFO("Longest SCC chain size is " << this->longestSccChainSize.value() << ".");
        if (needAdaptPrecision && &xLower == &xUpper) {
            // Each SCC requires (the same) increased precision eps' so that the overall accumulated error is longestChainSize * eps'
            setPrecisionDivisor(this->longestSccChainSize.value());
        }
    }

    // Handle the case where there is just one large SCC
    bool returnValue = true;
    if (this->sortedSccDecomposition->size() == 1) {
        if (auto const& scc = *this->sortedSccDecomposition->begin(); scc.size() == 1) {
            // Catch the trivial case where the whole system is just a single state.
            returnValue = solveTrivialScc(*scc.begin(), xUpper, bUpper);
            if (&xLower != &xUpper) {
                returnValue = solveTrivialScc(*scc.begin(), xLower, bLower) && returnValue;
            }
        } else {
            returnValue = solveFullyConnectedEquationSystem(sccSolverEnvironment, xLower, xUpper, bLower, bUpper);
        }
    } else {
        // Solve each SCC individually
        std::optional<storm::storage::BitVector> newRelevantValues;
        if (env.solver().topological().isExtendRelevantValues() && this->hasRelevantValues() &&
            this->sortedSccDecomposition->size() < this->A->getRowGroupCount()) {
            newRelevantValues = this->getRelevantValues();
            // Extend the relevant values towards those that have an incoming transition from another SCC
            std::vector<uint64_t> rowGroupToScc = this->sortedSccDecomposition->computeStateToSccIndexMap(this->A->getRowGroupCount());
            for (uint64_t rowGroup = 0; rowGroup < this->A->getRowGroupCount(); ++rowGroup) {
                auto currScc = rowGroupToScc[rowGroup];
                for (auto const& successor : this->A->getRowGroup(rowGroup)) {
                    if (rowGroupToScc[successor.getColumn()] != currScc) {
                        newRelevantValues->set(successor.getColumn(), true);
                    }
                }
            }
        }
        storm::storage::BitVector sccAsBitVector(xLower.size(), false);
        uint64_t sccIndex = 0;
        storm::utility::ProgressMeasurement progress("states");
        progress.setMaxCount(xLower.size());
        progress.startNewMeasurement(0);
        for (auto const& scc : *this->sortedSccDecomposition) {
            if (scc.size() == 1) {
                returnValue = solveTrivialScc(*scc.begin(), xUpper, bUpper);
                if (&xLower != &xUpper) {
                    returnValue = solveTrivialScc(*scc.begin(), xLower, bLower) && returnValue;
                }
            } else {
                sccAsBitVector.clear();
                for (auto const& state : scc) {
                    sccAsBitVector.set(state, true);
                }
                if (needAdaptPrecision && &xLower != &xUpper) {
                    // SCC require increased precision eps' based on their depth
                    setPrecisionDivisor(this->sortedSccDecomposition->getSccDepth(sccIndex) + 1);
                }
                returnValue = solveScc(sccSolverEnvironment, sccAsBitVector, xLower, xUpper, bLower, bUpper, newRelevantValues) && returnValue;
            }
            ++sccIndex;
            progress.updateProgress(sccIndex);
            if (storm::utility::resources::isTerminate()) {
                STORM_LOG_WARN("Topological solver aborted after analyzing " << sccIndex << "/" << this->sortedSccDecomposition->size() << " SCCs.");
                break;
            }
        }
    }

    if (!this->isCachingEnabled()) {
        clearCache();
    }

    return returnValue;
}

template<typename ValueType>
void TopologicalLinearEquationSolver<ValueType>::createSortedSccDecomposition(bool needSccDepths) const {
    // Obtain the scc decomposition
    this->sortedSccDecomposition = std::make_unique<storm::storage::StronglyConnectedComponentDecomposition<ValueType>>(
        *this->A, storm::storage::StronglyConnectedComponentDecompositionOptions().forceTopologicalSort().computeSccDepths(needSccDepths));
    if (needSccDepths) {
        this->longestSccChainSize = this->sortedSccDecomposition->getMaxSccDepth() + 1;
    }
}

template<typename ValueType>
bool TopologicalLinearEquationSolver<ValueType>::solveTrivialScc(uint64_t const& sccState, std::vector<ValueType>& globalX,
                                                                 std::vector<ValueType> const& globalB) const {
    ValueType& xi = globalX[sccState];
    xi = globalB[sccState];
    bool hasDiagonalEntry = false;
    ValueType denominator;
    for (auto const& entry : this->A->getRow(sccState)) {
        if (entry.getColumn() == sccState) {
            STORM_LOG_ASSERT(!storm::utility::isOne(entry.getValue()), "Diagonal entry of fix point system has value 1.");
            hasDiagonalEntry = true;
            denominator = storm::utility::one<ValueType>() - entry.getValue();
        } else {
            xi += entry.getValue() * globalX[entry.getColumn()];
        }
    }

    if (hasDiagonalEntry) {
        STORM_LOG_WARN_COND_DEBUG(
            storm::NumberTraits<ValueType>::IsExact || !storm::utility::isAlmostZero(denominator) || storm::utility::isZero(denominator),
            "State " << sccState << " has a selfloop with probability '1-(" << denominator << ")'. This could be an indication for numerical issues.");
        if (storm::utility::isZero(denominator)) {
            STORM_LOG_THROW(storm::utility::isZero(xi), storm::exceptions::InvalidOperationException, "The equation system has no solution.");
        } else {
            xi /= denominator;
        }
    }
    return true;
}

template<typename ValueType>
bool TopologicalLinearEquationSolver<ValueType>::solveFullyConnectedEquationSystem(storm::Environment const& sccSolverEnvironment,
                                                                                   std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                                   std::vector<ValueType> const& bLower,
                                                                                   std::vector<ValueType> const& bUpper) const {
    if (!this->sccSolver) {
        this->sccSolver = GeneralLinearEquationSolverFactory<ValueType>().create(sccSolverEnvironment);
        this->sccSolver->setCachingEnabled(true);
    }
    if (this->hasRelevantValues()) {
        this->sccSolver->setRelevantValues(this->getRelevantValues());
    }
    this->sccSolver->setBoundsFromOtherSolver(*this);
    if (this->sccSolver->getEquationProblemFormat(sccSolverEnvironment) == LinearEquationSolverProblemFormat::EquationSystem) {
        // Convert the matrix to an equation system. Note that we need to insert diagonal entries.
        storm::storage::SparseMatrix<ValueType> eqSysA(*this->A, true);
        eqSysA.convertToEquationSystem();
        this->sccSolver->setMatrix(std::move(eqSysA));
    } else {
        this->sccSolver->setMatrix(*this->A);
    }

    if (&xLower == &xUpper) {
        STORM_LOG_ASSERT(&bLower == &bUpper, "Provided different b-vectors but same x-vector.");
        return this->sccSolver->solveEquations(sccSolverEnvironment, xLower, bLower);
    } else {
        return this->sccSolver->solveEquationsSound(sccSolverEnvironment, xLower, xUpper, bLower, bUpper);
    }
}

template<typename ValueType>
bool TopologicalLinearEquationSolver<ValueType>::solveScc(storm::Environment const& sccSolverEnvironment, storm::storage::BitVector const& scc,
                                                          std::vector<ValueType>& xLowerGlobal, std::vector<ValueType>& xUpperGlobal,
                                                          std::vector<ValueType> const& bLowerGlobal, std::vector<ValueType> const& bUpperGlobal,
                                                          std::optional<storm::storage::BitVector> const& globalRelevantValues) const {
    // Set up the SCC solver
    if (!this->sccSolver) {
        this->sccSolver = GeneralLinearEquationSolverFactory<ValueType>().create(sccSolverEnvironment);
        this->sccSolver->setCachingEnabled(true);
    }
    if (globalRelevantValues) {
        this->sccSolver->setRelevantValues((*globalRelevantValues) % scc);
    }

    // Matrix
    bool asEquationSystem = this->sccSolver->getEquationProblemFormat(sccSolverEnvironment) == LinearEquationSolverProblemFormat::EquationSystem;
    storm::storage::SparseMatrix<ValueType> sccA = this->A->getSubmatrix(true, scc, scc, asEquationSystem);
    if (asEquationSystem) {
        sccA.convertToEquationSystem();
    }
    this->sccSolver->setMatrix(std::move(sccA));

    // x and b vectors
    auto prepareX = [&scc](std::vector<ValueType> const& globalX) { return storm::utility::vector::filterVector(globalX, scc); };
    auto prepareB = [&scc, this](std::vector<ValueType> const& globalX, std::vector<ValueType> const& globalB) {
        std::vector<ValueType> b;
        b.reserve(scc.getNumberOfSetBits());
        for (auto row : scc) {
            ValueType bi = globalB[row];
            for (auto const& entry : this->A->getRow(row)) {
                if (!scc.get(entry.getColumn())) {
                    bi += entry.getValue() * globalX[entry.getColumn()];
                }
            }
            b.push_back(std::move(bi));
        }
        return b;
    };

    // lower/upper bounds
    if (this->hasLowerBound(storm::solver::AbstractEquationSolver<ValueType>::BoundType::Global)) {
        this->sccSolver->setLowerBound(this->getLowerBound());
    } else if (this->hasLowerBound(storm::solver::AbstractEquationSolver<ValueType>::BoundType::Local)) {
        this->sccSolver->setLowerBounds(storm::utility::vector::filterVector(this->getLowerBounds(), scc));
    }
    if (this->hasUpperBound(storm::solver::AbstractEquationSolver<ValueType>::BoundType::Global)) {
        this->sccSolver->setUpperBound(this->getUpperBound());
    } else if (this->hasUpperBound(storm::solver::AbstractEquationSolver<ValueType>::BoundType::Local)) {
        this->sccSolver->setUpperBounds(storm::utility::vector::filterVector(this->getUpperBounds(), scc));
    }

    // std::cout << "rhs is " << storm::utility::vector::toString(sccB) << '\n';
    // std::cout << "x is " << storm::utility::vector::toString(sccX) << '\n';

    // Invoke scc solver

    bool returnvalue = false;
    if (&xLowerGlobal == &xUpperGlobal) {
        auto sccX = prepareX(xLowerGlobal);
        auto sccB = prepareB(xLowerGlobal, bLowerGlobal);
        returnvalue = this->sccSolver->solveEquations(sccSolverEnvironment, sccX, sccB);
        storm::utility::vector::setVectorValues(xLowerGlobal, scc, sccX);
    } else {
        auto sccXLower = prepareX(xLowerGlobal);
        auto sccXUpper = prepareX(xUpperGlobal);
        auto sccBLower = prepareB(xLowerGlobal, bLowerGlobal);
        auto sccBUpper = prepareB(xUpperGlobal, bUpperGlobal);
        returnvalue = this->sccSolver->solveEquationsSound(sccSolverEnvironment, sccXLower, sccXUpper, sccBLower, sccBUpper);
        storm::utility::vector::setVectorValues(xLowerGlobal, scc, sccXLower);
        storm::utility::vector::setVectorValues(xUpperGlobal, scc, sccXUpper);
    }
    return returnvalue;
}

template<typename ValueType>
LinearEquationSolverProblemFormat TopologicalLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const& env) const {
    return LinearEquationSolverProblemFormat::FixedPointSystem;
}

template<typename ValueType>
LinearEquationSolverRequirements TopologicalLinearEquationSolver<ValueType>::getRequirements(Environment const& env) const {
    // Return the requirements of the underlying solver
    return GeneralLinearEquationSolverFactory<ValueType>().getRequirements(getEnvironmentForUnderlyingSolver(env));
}

template<typename ValueType>
void TopologicalLinearEquationSolver<ValueType>::clearCache() const {
    sortedSccDecomposition.reset();
    longestSccChainSize = std::nullopt;
    sccSolver.reset();
    LinearEquationSolver<ValueType>::clearCache();
}

template<typename ValueType>
uint64_t TopologicalLinearEquationSolver<ValueType>::getMatrixRowCount() const {
    return this->A->getRowCount();
}

template<typename ValueType>
uint64_t TopologicalLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
    return this->A->getColumnCount();
}

template<typename ValueType>
std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> TopologicalLinearEquationSolverFactory<ValueType>::create(Environment const& env) const {
    return std::make_unique<storm::solver::TopologicalLinearEquationSolver<ValueType>>();
}

template<typename ValueType>
std::unique_ptr<LinearEquationSolverFactory<ValueType>> TopologicalLinearEquationSolverFactory<ValueType>::clone() const {
    return std::make_unique<TopologicalLinearEquationSolverFactory<ValueType>>(*this);
}

// Explicitly instantiate the linear equation solver.
template class TopologicalLinearEquationSolver<double>;
template class TopologicalLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
template class TopologicalLinearEquationSolver<storm::RationalNumber>;
template class TopologicalLinearEquationSolverFactory<storm::RationalNumber>;

template class TopologicalLinearEquationSolver<storm::RationalFunction>;
template class TopologicalLinearEquationSolverFactory<storm::RationalFunction>;

#endif
}  // namespace solver
}  // namespace storm

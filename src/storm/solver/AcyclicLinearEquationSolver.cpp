#include "storm/solver/AcyclicLinearEquationSolver.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/solver/helper/AcyclicSolverHelper.h"

#include "storm/utility/vector.h"

namespace storm {
namespace solver {

template<typename ValueType>
AcyclicLinearEquationSolver<ValueType>::AcyclicLinearEquationSolver() {
    // Intentionally left empty.
}

template<typename ValueType>
AcyclicLinearEquationSolver<ValueType>::AcyclicLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) {
    this->setMatrix(A);
}

template<typename ValueType>
AcyclicLinearEquationSolver<ValueType>::AcyclicLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) {
    this->setMatrix(std::move(A));
}

template<typename ValueType>
void AcyclicLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
    localA.reset();
    this->A = &A;
    clearCache();
}

template<typename ValueType>
void AcyclicLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
    localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
    this->A = localA.get();
    clearCache();
}

template<typename ValueType>
uint64_t AcyclicLinearEquationSolver<ValueType>::getMatrixRowCount() const {
    return this->A->getRowCount();
}

template<typename ValueType>
uint64_t AcyclicLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
    return this->A->getColumnCount();
}

template<typename ValueType>
bool AcyclicLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                    std::vector<ValueType> const& bLower, std::vector<ValueType> const& bUpper) const {
    if (&xLower != &xUpper) {
        bool ret = internalSolveEquations(env, xLower, xLower, bLower, bLower);
        if (&bLower == &bUpper) {
            xUpper = xLower;  // same b vector, no need to recompute. Just copy over values
        } else {
            ret = internalSolveEquations(env, xUpper, xUpper, bUpper, bUpper) && ret;
        }
        return ret;
    }
    STORM_LOG_ASSERT(&bLower == &bUpper, "Invalid function call: different lower/upper b vector but same result vector.");

    auto& x = xLower;
    auto& b = bLower;

    STORM_LOG_ASSERT(x.size() == this->A->getRowGroupCount(), "Provided x-vector has invalid size.");
    STORM_LOG_ASSERT(b.size() == this->A->getRowCount(), "Provided b-vector has invalid size.");

    if (!multiplier) {
        // We have not allocated cache memory, yet
        rowOrdering = helper::computeTopologicalGroupOrdering(*this->A);
        if (!rowOrdering) {
            // It is not required to reorder the elements.
            this->multiplier = storm::solver::MultiplierFactory<ValueType>().create(env, *this->A);
        } else {
            bFactors.clear();
            orderedMatrix = helper::createReorderedMatrix(*this->A, *rowOrdering, bFactors);
            this->multiplier = storm::solver::MultiplierFactory<ValueType>().create(env, *orderedMatrix);
        }
        auxiliaryRowVector = std::vector<ValueType>(this->A->getRowCount());
        auxiliaryRowVector2 = std::vector<ValueType>(this->A->getRowCount());
    }

    std::vector<ValueType>* xPtr = &x;
    std::vector<ValueType> const* bPtr = &b;
    if (rowOrdering) {
        STORM_LOG_ASSERT(rowOrdering->size() == b.size(), "b-vector has unexpected size.");
        auxiliaryRowVector->resize(b.size());
        storm::utility::vector::selectVectorValues(*auxiliaryRowVector, *rowOrdering, b);
        for (auto const& bFactor : bFactors) {
            (*auxiliaryRowVector)[bFactor.first] *= bFactor.second;
        }
        bPtr = &auxiliaryRowVector.get();
        xPtr = &auxiliaryRowVector2.get();
    }

    this->multiplier->multiplyGaussSeidel(env, *xPtr, bPtr, true);

    if (rowOrdering) {
        for (uint64_t newRow = 0; newRow < x.size(); ++newRow) {
            x[(*rowOrdering)[newRow]] = (*xPtr)[newRow];
        }
    }

    if (!this->isCachingEnabled()) {
        this->clearCache();
    }
    return true;
}

template<typename ValueType>
LinearEquationSolverProblemFormat AcyclicLinearEquationSolver<ValueType>::getEquationProblemFormat(storm::Environment const& env) const {
    return LinearEquationSolverProblemFormat::FixedPointSystem;
}

template<typename ValueType>
LinearEquationSolverRequirements AcyclicLinearEquationSolver<ValueType>::getRequirements(Environment const& env) const {
    // Return the requirements of the underlying solver
    LinearEquationSolverRequirements requirements;
    requirements.requireAcyclic();
    return requirements;
}

template<typename ValueType>
bool AcyclicLinearEquationSolver<ValueType>::supportsSolveEquationsSoundBounds(Environment const&, bool) const {
    return true;
}

template<typename ValueType>
void AcyclicLinearEquationSolver<ValueType>::clearCache() const {
    multiplier.reset();
    orderedMatrix = boost::none;
    rowOrdering = boost::none;
    auxiliaryRowVector = boost::none;
    auxiliaryRowVector2 = boost::none;
    bFactors.clear();
}

// Explicitly instantiate the min max linear equation solver.
template class AcyclicLinearEquationSolver<double>;

#ifdef STORM_HAVE_CARL
template class AcyclicLinearEquationSolver<storm::RationalNumber>;
template class AcyclicLinearEquationSolver<storm::RationalFunction>;
#endif
}  // namespace solver
}  // namespace storm

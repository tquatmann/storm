#include "storm/solver/EigenLinearEquationSolver.h"

#include "storm/adapters/EigenAdapter.h"

#include "storm/adapters/RationalFunctionAdapter.h"

#include "storm/environment/solver/EigenSolverEnvironment.h"

#include "storm/exceptions/InvalidSettingsException.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

namespace storm {
namespace solver {

template<typename ValueType>
EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver() {
    // Intentionally left empty.
}

template<typename ValueType>
EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) {
    this->setMatrix(A);
}

template<typename ValueType>
EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) {
    this->setMatrix(std::move(A));
}

template<typename ValueType>
void EigenLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
    eigenA = storm::adapters::EigenAdapter::toEigenSparseMatrix<ValueType>(A);
    this->clearCache();
}

template<typename ValueType>
void EigenLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
    // Take ownership of the matrix so it is destroyed after we have translated it to Eigen's format.
    storm::storage::SparseMatrix<ValueType> localA(std::move(A));
    this->setMatrix(localA);
    this->clearCache();
}

template<typename ValueType>
EigenLinearEquationSolverMethod EigenLinearEquationSolver<ValueType>::getMethod(Environment const& env, bool isExactMode) const {
    // Adjust the method if none was specified and we are using rational numbers.
    auto method = env.solver().eigen().getMethod();

    if (isExactMode && method != EigenLinearEquationSolverMethod::SparseLU) {
        if (env.solver().eigen().isMethodSetFromDefault()) {
            STORM_LOG_INFO("Selecting 'SparseLU' as the solution technique to guarantee exact results.");
        } else {
            STORM_LOG_WARN("The selected solution method does not guarantee exact results. Falling back to SparseLU.");
        }
        method = EigenLinearEquationSolverMethod::SparseLU;
    } else if (env.solver().isForceSoundness() && method != EigenLinearEquationSolverMethod::SparseLU) {
        if (env.solver().eigen().isMethodSetFromDefault()) {
            STORM_LOG_INFO(
                "Selecting 'SparseLU' as the solution technique to guarantee sound results. If you want to override this, please explicitly specify a "
                "different method.");
            method = EigenLinearEquationSolverMethod::SparseLU;
        } else {
            STORM_LOG_WARN("The selected solution method does not guarantee sound results.");
        }
    }
    return method;
}

template<typename ValueType>
bool EigenLinearEquationSolver<ValueType>::supportsSolveEquationsSoundBounds(Environment const& env, bool) const {
    return getMethod(env, storm::NumberTraits<ValueType>::IsExact || env.solver().isForceExact()) == EigenLinearEquationSolverMethod::SparseLU;
}

template<typename ValueType>
bool EigenLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& xLower, std::vector<ValueType>& xUpper,
                                                                  std::vector<ValueType> const& bLower, std::vector<ValueType> const& bUpper) const {
    // Determine the method and whether it is exact (based on environment and ValueType)
    bool const isExact = storm::NumberTraits<ValueType>::IsExact || env.solver().isForceExact();
    auto solutionMethod = getMethod(env, isExact);

    // Map the input vectors to Eigen's format.
    // Only consider xLower/bLower for now. If necessary, upper bounds are computed in the SparseLU case below.
    auto eigenVec = [](auto& vec) { return Eigen::Matrix<ValueType, Eigen::Dynamic, 1>::Map(vec.data(), vec.size()); };
    auto eigenX = eigenVec(xLower);
    auto eigenB = eigenVec(bLower);

    if (solutionMethod == EigenLinearEquationSolverMethod::SparseLU) {
        STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with sparse LU factorization (Eigen library).");
        Eigen::SparseLU<Eigen::SparseMatrix<ValueType>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(*this->eigenA);
        solver._solve_impl(eigenB, eigenX);
        if (solver.info() == Eigen::ComputationInfo::Success && &xLower != &xUpper) {
            if (&bLower == &bUpper) {
                xUpper = xLower;  // same b vector, no need to recompute. Just copy over values
            } else {
                // Need to compute upper bounds as well.
                auto eigenXUpper = eigenVec(xUpper);
                auto eigenBUpper = eigenVec(bUpper);
                solver._solve_impl(eigenBUpper, eigenXUpper);
            }
        }
        return solver.info() == Eigen::ComputationInfo::Success;
    }
    if constexpr (storm::NumberTraits<ValueType>::IsExact) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidSettingsException, "Unable to perform numerical computations on exact value type.");
    } else {
        STORM_LOG_ASSERT(!isExact, "Unexpected method selected for exact computation.");
        STORM_LOG_THROW(&xLower == &xUpper, storm::exceptions::NotImplementedException,
                        "EigenLinearEquationSolver in non-exact mode does not support solving with sound lower and upper bounds.");
        STORM_LOG_ASSERT(&bLower == &bUpper, "Invalid function call: different lower/upper b vector but same result vector.");

        bool converged = false;
        uint64_t numberOfIterations = 0;
        Eigen::Index maxIter = std::numeric_limits<Eigen::Index>::max();
        if (env.solver().eigen().getMaximalNumberOfIterations() < static_cast<uint64_t>(maxIter)) {
            maxIter = env.solver().eigen().getMaximalNumberOfIterations();
        }
        uint64_t restartThreshold = env.solver().eigen().getRestartThreshold();
        ValueType precision = storm::utility::convertNumber<ValueType>(env.solver().eigen().getPrecision());
        EigenLinearEquationSolverPreconditioner preconditioner = env.solver().eigen().getPreconditioner();
        if (solutionMethod == EigenLinearEquationSolverMethod::Bicgstab) {
            if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with BiCGSTAB with Ilu preconditioner (Eigen library).");

                Eigen::BiCGSTAB<Eigen::SparseMatrix<ValueType>, Eigen::IncompleteLUT<ValueType>> solver;
                solver.compute(*this->eigenA);
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with BiCGSTAB with Diagonal preconditioner (Eigen library).");

                Eigen::BiCGSTAB<Eigen::SparseMatrix<ValueType>, Eigen::DiagonalPreconditioner<ValueType>> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with BiCGSTAB with identity preconditioner (Eigen library).");

                Eigen::BiCGSTAB<Eigen::SparseMatrix<ValueType>, Eigen::IdentityPreconditioner> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                numberOfIterations = solver.iterations();
                converged = solver.info() == Eigen::ComputationInfo::Success;
            }
        } else if (solutionMethod == EigenLinearEquationSolverMethod::DGmres) {
            if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with DGMRES with Ilu preconditioner (Eigen library).");
                Eigen::DGMRES<Eigen::SparseMatrix<ValueType>, Eigen::IncompleteLUT<ValueType>> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with DGMRES with Diagonal preconditioner (Eigen library).");

                Eigen::DGMRES<Eigen::SparseMatrix<ValueType>, Eigen::DiagonalPreconditioner<ValueType>> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with DGMRES with identity preconditioner (Eigen library).");

                Eigen::DGMRES<Eigen::SparseMatrix<ValueType>, Eigen::IdentityPreconditioner> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            }

        } else if (solutionMethod == EigenLinearEquationSolverMethod::Gmres) {
            if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with GMRES with Ilu preconditioner (Eigen library).");

                Eigen::GMRES<Eigen::SparseMatrix<ValueType>, Eigen::IncompleteLUT<ValueType>> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with GMRES with Diagonal preconditioner (Eigen library).");

                Eigen::GMRES<Eigen::SparseMatrix<ValueType>, Eigen::DiagonalPreconditioner<ValueType>> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            } else {
                STORM_LOG_INFO("Solving linear equation system (" << xLower.size() << " rows) with GMRES with identity preconditioner (Eigen library).");

                Eigen::GMRES<Eigen::SparseMatrix<ValueType>, Eigen::IdentityPreconditioner> solver;
                solver.setTolerance(precision);
                solver.setMaxIterations(maxIter);
                solver.set_restart(restartThreshold);
                solver.compute(*this->eigenA);
                eigenX = solver.solveWithGuess(eigenB, eigenX);
                converged = solver.info() == Eigen::ComputationInfo::Success;
                numberOfIterations = solver.iterations();
            }
        }

        // Make sure that all results conform to the (global) bounds.
        storm::utility::vector::clip(xLower, this->lowerBound, this->upperBound);

        // Check if the solver converged and issue a warning otherwise.
        if (converged) {
            STORM_LOG_INFO("Iterative solver converged after " << numberOfIterations << " iterations.");
            return true;
        } else {
            STORM_LOG_WARN("Iterative solver did not converge.");
            return false;
        }
    }
}

template<typename ValueType>
LinearEquationSolverProblemFormat EigenLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const&) const {
    return LinearEquationSolverProblemFormat::EquationSystem;
}

template<typename ValueType>
uint64_t EigenLinearEquationSolver<ValueType>::getMatrixRowCount() const {
    return eigenA->rows();
}

template<typename ValueType>
uint64_t EigenLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
    return eigenA->cols();
}

template<typename ValueType>
std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> EigenLinearEquationSolverFactory<ValueType>::create(Environment const&) const {
    return std::make_unique<storm::solver::EigenLinearEquationSolver<ValueType>>();
}

template<typename ValueType>
std::unique_ptr<LinearEquationSolverFactory<ValueType>> EigenLinearEquationSolverFactory<ValueType>::clone() const {
    return std::make_unique<EigenLinearEquationSolverFactory<ValueType>>(*this);
}

template class EigenLinearEquationSolver<double>;
template class EigenLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
template class EigenLinearEquationSolver<storm::RationalNumber>;
template class EigenLinearEquationSolver<storm::RationalFunction>;

template class EigenLinearEquationSolverFactory<storm::RationalNumber>;
template class EigenLinearEquationSolverFactory<storm::RationalFunction>;
#endif
}  // namespace solver
}  // namespace storm

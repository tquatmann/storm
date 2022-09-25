#include "storm/solver/LpMinMaxLinearEquationSolver.h"

#include "storm/environment/solver/MinMaxLpSolverEnvironment.h"
#include "storm/environment/solver/MinMaxSolverEnvironment.h"
#include "storm/exceptions/InvalidEnvironmentException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/storage/expressions/RationalLiteralExpression.h"
#include "storm/storage/expressions/VariableExpression.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

// Define to write LP to a file
#undef STORM_LPMINMAX_DEBUG

namespace storm {
namespace solver {

template<typename ValueType>
LpMinMaxLinearEquationSolver<ValueType>::LpMinMaxLinearEquationSolver(std::unique_ptr<storm::utility::solver::LpSolverFactory<ValueType>>&& lpSolverFactory)
    : lpSolverFactory(std::move(lpSolverFactory)) {
    // Intentionally left empty.
}

template<typename ValueType>
LpMinMaxLinearEquationSolver<ValueType>::LpMinMaxLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A,
                                                                      std::unique_ptr<storm::utility::solver::LpSolverFactory<ValueType>>&& lpSolverFactory)
    : StandardMinMaxLinearEquationSolver<ValueType>(A), lpSolverFactory(std::move(lpSolverFactory)) {
    // Intentionally left empty.
}

template<typename ValueType>
LpMinMaxLinearEquationSolver<ValueType>::LpMinMaxLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A,
                                                                      std::unique_ptr<storm::utility::solver::LpSolverFactory<ValueType>>&& lpSolverFactory)
    : StandardMinMaxLinearEquationSolver<ValueType>(std::move(A)), lpSolverFactory(std::move(lpSolverFactory)) {
    // Intentionally left empty.
}

template<typename ValueType>
bool LpMinMaxLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, OptimizationDirection dir, std::vector<ValueType>& x,
                                                                     std::vector<ValueType> const& b) const {
    STORM_LOG_THROW(env.solver().minMax().getMethod() == MinMaxMethod::LinearProgramming, storm::exceptions::InvalidEnvironmentException,
                    "This min max solver does not support the selected technique.");
    bool const doOptimize = env.solver().minMax().getMinMaxLpSolverEnvironment().getOptimizeFeasibilityQueries() || !this->hasGlobalBound();
    // Set up the LP solver
    std::unique_ptr<storm::solver::LpSolver<ValueType>> solver = lpSolverFactory->create("");
    solver->setOptimizationDirection(invert(dir));
    bool const optimizeOnlyRelevant = this->hasRelevantValues() && env.solver().minMax().getMinMaxLpSolverEnvironment().getOptimizeOnlyForInitialState();
    STORM_LOG_DEBUG("Optimize only for relevant state requested:" << env.solver().minMax().getMinMaxLpSolverEnvironment().getOptimizeOnlyForInitialState());
    if (optimizeOnlyRelevant) {
        STORM_LOG_TRACE("Relevant values " << this->getRelevantValues());
    } else if (!this->hasRelevantValues()) {
        STORM_LOG_DEBUG("No relevant values set! Optimizing over all states.");
    }
    bool const useBounds = env.solver().minMax().getMinMaxLpSolverEnvironment().getUseNonTrivialBounds();
    // Create a variable for each row group
    std::vector<storm::expressions::Expression> variableExpressions;
    variableExpressions.reserve(this->A->getRowGroupCount());
    for (uint64_t rowGroup = 0; rowGroup < this->A->getRowGroupCount(); ++rowGroup) {
        ValueType objValue = !doOptimize || (optimizeOnlyRelevant && !this->getRelevantValues().get(rowGroup)) ? storm::utility::zero<ValueType>()
                                                                                                               : storm::utility::one<ValueType>();
        if (useBounds && this->hasLowerBound()) {
            ValueType lowerBound = this->getLowerBound(rowGroup);
            if (this->hasUpperBound()) {
                ValueType upperBound = this->getUpperBound(rowGroup);
                if (lowerBound == upperBound) {
                    // Some solvers (like glpk) don't support variables with bounds [x,x]. We therefore just use a constant instead. This should be more
                    // efficient anyways.
                    variableExpressions.push_back(solver->getConstant(lowerBound));
                } else {
                    STORM_LOG_ASSERT(lowerBound <= upperBound,
                                     "Lower Bound at row group " << rowGroup << " is " << lowerBound << " which exceeds the upper bound " << upperBound << ".");
                    variableExpressions.emplace_back(solver->addBoundedContinuousVariable("x" + std::to_string(rowGroup), lowerBound, upperBound, objValue));
                }
            } else {
                variableExpressions.emplace_back(solver->addLowerBoundedContinuousVariable("x" + std::to_string(rowGroup), lowerBound, objValue));
            }
        } else {
            if (useBounds && this->hasUpperBound()) {
                variableExpressions.emplace_back(
                    solver->addUpperBoundedContinuousVariable("x" + std::to_string(rowGroup), this->getUpperBound(rowGroup), objValue));
            } else {
                variableExpressions.emplace_back(solver->addUnboundedContinuousVariable("x" + std::to_string(rowGroup), objValue));
            }
        }
    }
    solver->update();
    STORM_LOG_DEBUG("Use eq if there is a single action: " << env.solver().minMax().getMinMaxLpSolverEnvironment().getUseEqualityForSingleActions());
    bool const useEqualityForSingleAction = env.solver().minMax().getMinMaxLpSolverEnvironment().getUseEqualityForSingleActions();

    if (this->hasGlobalBound()) {
        STORM_LOG_DEBUG("Global bound set. Add constraints for relevant values...");
        storm::expressions::Expression thresholdConstant = solver->getConstant(this->globalBound->constraintValue);
        storm::expressions::RelationType reltype;
        if (this->globalBound->lowerBound) {
            if (this->globalBound->strict) {
                reltype = storm::expressions::RelationType::Greater;
            } else {
                reltype = storm::expressions::RelationType::GreaterOrEqual;
            }
        } else {
            if (this->globalBound->strict) {
                reltype = storm::expressions::RelationType::Less;
            } else {
                reltype = storm::expressions::RelationType::LessOrEqual;
            }
        }
        for (uint64_t entry : this->getRelevantValues()) {
            storm::expressions::Expression constraint =
                storm::expressions::makeBinaryRelationExpression(variableExpressions[entry], thresholdConstant, reltype);
            STORM_LOG_TRACE("Global bound added: " << constraint);
            solver->addConstraint("", constraint);
        }
    }

    // Add a constraint for each row
    for (uint64_t rowGroup = 0; rowGroup < this->A->getRowGroupCount(); ++rowGroup) {
        // The rowgroup refers to the state number
        uint64_t rowIndex, rowGroupEnd;
        if (this->choiceFixedForRowGroup && this->choiceFixedForRowGroup.get()[rowGroup]) {
            rowIndex = this->A->getRowGroupIndices()[rowGroup] + this->getInitialScheduler()[rowGroup];
            rowGroupEnd = rowIndex + 1;
        } else {
            rowIndex = this->A->getRowGroupIndices()[rowGroup];
            rowGroupEnd = this->A->getRowGroupIndices()[rowGroup + 1];
        }
        bool const singleAction = (rowIndex + 1 == rowGroupEnd);
        bool const useEquality = useEqualityForSingleAction && singleAction;
        for (; rowIndex < rowGroupEnd; ++rowIndex) {
            auto row = this->A->getRow(rowIndex);
            std::vector<storm::expressions::Expression> summands;
            summands.reserve(1 + row.getNumberOfEntries());
            summands.push_back(solver->getConstant(b[rowIndex]));
            for (auto const& entry : row) {
                summands.push_back(solver->getConstant(entry.getValue()) * variableExpressions[entry.getColumn()]);
            }
            storm::expressions::Expression rowConstraint = storm::expressions::sum(summands);
            if (useEquality) {
                rowConstraint = variableExpressions[rowGroup] == rowConstraint;
            } else if (minimize(dir)) {
                rowConstraint = variableExpressions[rowGroup] <= rowConstraint;
            } else {
                rowConstraint = variableExpressions[rowGroup] >= rowConstraint;
            }
            solver->addConstraint("", rowConstraint);
        }
    }

#ifdef STORM_LPMINMAX_DEBUG
    solver->writeModelToFile("debug.lp");
#endif

    // Invoke optimization
    STORM_LOG_TRACE("Run solver...");
    solver->optimize();
    STORM_LOG_TRACE("...done.");

    bool infeasible = solver->isInfeasible();
    if (!this->hasGlobalBound()) {
        STORM_LOG_THROW(!infeasible, storm::exceptions::UnexpectedException, "The MinMax equation system is infeasible.");
    } else {
        // Explicitly nothing to be done here.
    }
    if (!infeasible) {
        STORM_LOG_THROW(!solver->isUnbounded(), storm::exceptions::UnexpectedException, "The MinMax equation system is unbounded.");
        STORM_LOG_THROW(solver->isOptimal(), storm::exceptions::UnexpectedException, "Unable to find optimal solution for MinMax equation system.");
    }

    if (!infeasible) {
        // write the solution into the solution vector
        STORM_LOG_ASSERT(x.size() == variableExpressions.size(), "Dimension of x-vector does not match number of varibales.");
        auto xIt = x.begin();
        auto vIt = variableExpressions.begin();
        for (; xIt != x.end(); ++xIt, ++vIt) {
            auto const& vBaseExpr = vIt->getBaseExpression();
            if (vBaseExpr.isVariable()) {
                *xIt = solver->getContinuousValue(vBaseExpr.asVariableExpression().getVariable());
            } else {
                STORM_LOG_ASSERT(vBaseExpr.isRationalLiteralExpression(), "Variable expression has unexpected type.");
                *xIt = storm::utility::convertNumber<ValueType>(vBaseExpr.asRationalLiteralExpression().getValue());
            }
        }

        // If requested, we store the scheduler for retrieval.
        if (this->isTrackSchedulerSet()) {
            this->schedulerChoices = std::vector<uint_fast64_t>(this->A->getRowGroupCount());
            for (uint64_t rowGroup = 0; rowGroup < this->A->getRowGroupCount(); ++rowGroup) {
                if (!this->choiceFixedForRowGroup || !this->choiceFixedForRowGroup.get()[rowGroup]) {
                    // Only update scheduler choice for the states that don't have a fixed choice
                    uint64_t row = this->A->getRowGroupIndices()[rowGroup];
                    uint64_t optimalChoiceIndex = 0;
                    uint64_t currChoice = 0;
                    ValueType optimalGroupValue = this->A->multiplyRowWithVector(row, x) + b[row];
                    for (++row, ++currChoice; row < this->A->getRowGroupIndices()[rowGroup + 1]; ++row, ++currChoice) {
                        ValueType rowValue = this->A->multiplyRowWithVector(row, x) + b[row];
                        if ((minimize(dir) && rowValue < optimalGroupValue) || (maximize(dir) && rowValue > optimalGroupValue)) {
                            optimalGroupValue = rowValue;
                            optimalChoiceIndex = currChoice;
                        }
                    }
                    this->schedulerChoices.get()[rowGroup] = optimalChoiceIndex;
                }
            }
        }
    } else {
        STORM_LOG_ASSERT(this->hasGlobalBound(), "Infeasbile solutions only occur due to a global bound.");
        // result is infeasible.
        // we report a result that violates the bound.

        ValueType res = this->globalBound->constraintValue;
        if (!this->globalBound->strict) {
            if (!this->globalBound->lowerBound) {
                res -= 1;
            } else {
                res += 1;
            }
        }
        for (auto xIt = x.begin(); xIt != x.end(); ++xIt) {
            *xIt = res;
        }

        // if tracking a scheduler was set, we report that we cannot find a scheduler as the bound is violated.
        if (this->isTrackSchedulerSet()) {
            STORM_LOG_WARN("Cannot track a scheduler if the bound is violated.");
        }
    }

    // Why always return true?
    return true;
}

template<typename ValueType>
void LpMinMaxLinearEquationSolver<ValueType>::clearCache() const {
    StandardMinMaxLinearEquationSolver<ValueType>::clearCache();
}

template<typename ValueType>
MinMaxLinearEquationSolverRequirements LpMinMaxLinearEquationSolver<ValueType>::getRequirements(
    Environment const& env, boost::optional<storm::solver::OptimizationDirection> const& direction, bool const& hasInitialScheduler) const {
    MinMaxLinearEquationSolverRequirements requirements;

    // In case we need to retrieve a scheduler, the solution has to be unique
    if (!this->hasUniqueSolution() && (env.solver().minMax().isForceRequireUnique() || this->isTrackSchedulerSet())) {
        requirements.requireUniqueSolution();
    }

    if (env.solver().minMax().getMinMaxLpSolverEnvironment().getUseNonTrivialBounds()) {
        requirements.requireBounds(false);
    }

    return requirements;
}

template class LpMinMaxLinearEquationSolver<double>;

#ifdef STORM_HAVE_CARL
template class LpMinMaxLinearEquationSolver<storm::RationalNumber>;
#endif
}  // namespace solver
}  // namespace storm

#include <algorithm>
#include <iterator>

#include "storm/modelchecker/helper/conditional/ConditionalHelper.h"

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/environment/modelchecker/ModelCheckerEnvironment.h"
#include "storm/environment/solver/MinMaxSolverEnvironment.h"
#include "storm/modelchecker/prctl/helper/SparseMdpPrctlHelper.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/solver/SolveGoal.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/transformer/EndComponentEliminator.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/KwekMehlhorn.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/constants.h"
#include "storm/utility/graph.h"
#include "storm/utility/logging.h"
#include "storm/utility/macros.h"

#include "storm/modelchecker/prctl/helper/MDPModelCheckingHelperReturnType.h"

#include "storm/exceptions/NotImplementedException.h"
#include "storm/exceptions/NotSupportedException.h"

namespace storm::modelchecker {

namespace internal {

template<typename ValueType>
boost::optional<typename storm::transformer::EndComponentEliminator<ValueType>::EndComponentEliminatorReturnType> eliminateEndComponents(
    storm::storage::BitVector possibleEcStates, bool addRowAtRepresentativeState, std::optional<uint64_t> representativeRowEntry,
    storm::storage::SparseMatrix<ValueType>& matrix, uint64_t& initialState, storm::storage::BitVector& rowsWithSum1, std::vector<ValueType>& rowValues1,
    storm::OptionalRef<std::vector<ValueType>> rowValues2 = {}) {
    storm::storage::MaximalEndComponentDecomposition<ValueType> ecs(matrix, matrix.transpose(true), possibleEcStates, rowsWithSum1);
    if (ecs.empty()) {
        return boost::none;  // nothing to do
    }

    storm::storage::BitVector allRowGroups(matrix.getRowGroupCount(), true);
    auto ecElimResult = storm::transformer::EndComponentEliminator<ValueType>::transform(
        matrix, ecs, allRowGroups, addRowAtRepresentativeState ? allRowGroups : ~allRowGroups, representativeRowEntry.has_value());

    // Update matrix
    matrix = std::move(ecElimResult.matrix);
    if (addRowAtRepresentativeState && representativeRowEntry) {
        auto const columnIndex = ecElimResult.oldToNewStateMapping[*representativeRowEntry];
        for (auto representativeRowIndex : ecElimResult.sinkRows) {
            auto row = matrix.getRow(representativeRowIndex);
            STORM_LOG_ASSERT(row.getNumberOfEntries() == 1, "unexpected number of entries in representative row.");
            auto& entry = *row.begin();
            entry.setColumn(columnIndex);
        }
    }

    // update vectors
    auto updateRowValue = [&ecElimResult](std::vector<ValueType>& rowValues) {
        std::vector<ValueType> newRowValues;
        newRowValues.reserve(ecElimResult.newToOldRowMapping.size());
        for (auto oldRowIndex : ecElimResult.newToOldRowMapping) {
            newRowValues.push_back(rowValues[oldRowIndex]);
        }
        rowValues = std::move(newRowValues);
        STORM_LOG_ASSERT(
            std::all_of(ecElimResult.sinkRows.begin(), ecElimResult.sinkRows.end(), [&rowValues](auto i) { return storm::utility::isZero(rowValues[i]); }),
            "Sink rows are expected to have zero value");
    };
    updateRowValue(rowValues1);
    if (rowValues2) {
        updateRowValue(*rowValues2);
    }

    // update initial state
    initialState = ecElimResult.oldToNewStateMapping[initialState];

    // update bitvector
    storm::storage::BitVector newRowsWithSum1(ecElimResult.newToOldRowMapping.size(), true);
    uint64_t newRowIndex = 0;
    for (auto oldRowIndex : ecElimResult.newToOldRowMapping) {
        if ((addRowAtRepresentativeState && !representativeRowEntry.has_value() && ecElimResult.sinkRows.get(newRowIndex)) || !rowsWithSum1.get(oldRowIndex)) {
            newRowsWithSum1.set(newRowIndex, false);
        }
        ++newRowIndex;
    }
    rowsWithSum1 = std::move(newRowsWithSum1);

    return ecElimResult;
}

template<typename ValueType>
struct SolverResult {
    SolverResult(ValueType initialStateValue) : initialStateValue(initialStateValue) {
        // Intentionally left empty.
    }

    bool hasScheduler() const {
        return static_cast<bool>(scheduler);
    }

    ValueType initialStateValue;
    boost::optional<std::vector<uint64_t>> scheduler;
};

template<typename ValueType, typename SolutionType = ValueType>
typename internal::SolverResult<ValueType> solveMinMaxEquationSystem(storm::Environment const& env, storm::storage::SparseMatrix<ValueType> const& matrix,
                                                                     std::vector<ValueType> const& rowValues, storm::storage::BitVector const& rowsWithSum1,
                                                                     storm::solver::OptimizationDirection const dir, uint64_t const initialState) {
    // Initialize the solution vector.
    std::vector<SolutionType> x(matrix.getRowGroupCount(), storm::utility::zero<ValueType>());

    // Set up the solver.
    auto solver = storm::solver::GeneralMinMaxLinearEquationSolverFactory<ValueType, SolutionType>().create(env, matrix);
    solver->setOptimizationDirection(dir);
    solver->setRequirementsChecked();
    solver->setHasUniqueSolution(true);
    solver->setHasNoEndComponents(true);
    solver->setLowerBound(storm::utility::zero<ValueType>());
    solver->setUpperBound(storm::utility::one<ValueType>());
    solver->setTrackScheduler(true);

    // Solve the corresponding system of equations.
    solver->solveEquations(env, x, rowValues);

    SolverResult<SolutionType> result(x[initialState]);
    result.scheduler = std::move(solver->getSchedulerChoices());

    return result;
}

/*!
 * Computes the reachability probabilities for the given target states and inserts all non-zero values into the given map.
 * @note This code is optimized for cases where not all states are reachable from the initial states.
 */
template<typename ValueType>
void computeReachabilityProbabilities(Environment const& env, std::map<uint64_t, ValueType>& nonZeroResults, storm::solver::OptimizationDirection const dir,
                                      storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::BitVector const& initialStates,
                                      storm::storage::BitVector const& allowedStates, storm::storage::BitVector const& targetStates) {
    if (initialStates.empty()) {  // nothing to do
        return;
    }
    auto const reachableStates = storm::utility::graph::getReachableStates(transitionMatrix, initialStates, allowedStates, targetStates);
    auto const subTargets = targetStates % reachableStates;
    // Catch the case where no target is reachable from an initial state. In this case, there is nothing to do since all probabilities are zero.
    if (subTargets.empty()) {
        return;
    }
    auto const subInits = initialStates % reachableStates;
    auto const submatrix = transitionMatrix.getSubmatrix(true, reachableStates, reachableStates);
    auto const subResult = helper::SparseMdpPrctlHelper<ValueType, ValueType>::computeUntilProbabilities(
        env, storm::solver::SolveGoal<ValueType>(dir, subInits), submatrix, submatrix.transpose(true), storm::storage::BitVector(subTargets.size(), true),
        subTargets, false, true);
    auto origInitIt = initialStates.begin();
    for (auto subInit : subInits) {
        auto const& val = subResult.values[subInit];
        if (!storm::utility::isZero(val)) {
            nonZeroResults.emplace(*origInitIt, val);
        }
        ++origInitIt;
    }
}

/*!
 * Computes the reachability probabilities for the given target states and inserts all non-zero values into the given map.
 * TODO we can optimize this by only considering reachable states as well but I don't have the strength to do it now
 */
template<typename ValueType>
std::unique_ptr<storm::storage::Scheduler<ValueType>> computeReachabilityProbabilitiesAndScheduler(
    Environment const& env, std::map<uint64_t, ValueType>& nonZeroResults, storm::solver::OptimizationDirection const dir,
    storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::BitVector const& initialStates,
    storm::storage::BitVector const& allowedStates, storm::storage::BitVector const& targetStates) {
    if (initialStates.empty()) {  // nothing to do
        return nullptr;
    }

    auto result = helper::SparseMdpPrctlHelper<ValueType, ValueType>::computeUntilProbabilities(
        env, storm::solver::SolveGoal<ValueType>(dir, initialStates), transitionMatrix, transitionMatrix.transpose(true),
        storm::storage::BitVector(targetStates.size(), true), targetStates, false, true);
    for (auto initState : initialStates) {
        auto const& val = result.values[initState];
        if (!storm::utility::isZero(val)) {
            nonZeroResults.emplace(initState, val);
        }
    }

    return std::move(result.scheduler);
}

/*!
 * Uses the precomputed the reachability probabilities for the given target states and inserts all non-zero values into the given map.
 * @note This code is optimized for cases where not all states are reachable from the initial states.
 */
template<typename ValueType>
void computeReachabilityProbabilities(Environment const& env, std::map<uint64_t, ValueType>& nonZeroResults, storm::solver::OptimizationDirection const dir,
                                      storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::BitVector const& initialStates,
                                      storm::storage::BitVector const& allowedStates, storm::storage::BitVector const& targetStates,
                                      std::vector<double> const& precomputedReachProbs) {
    if (initialStates.empty()) {  // nothing to do
        return;
    }
    auto const reachableStates = storm::utility::graph::getReachableStates(transitionMatrix, initialStates, allowedStates, targetStates);
    auto const subTargets = targetStates % reachableStates;
    // Catch the case where no target is reachable from an initial state. In this case, there is nothing to do since all probabilities are zero.
    if (subTargets.empty()) {
        return;
    }

    for (auto initState : initialStates) {
        auto const& val = precomputedReachProbs[initState];
        if (!storm::utility::isZero(val)) {
            nonZeroResults.emplace(initState, val);
        }
    }
}

template<typename ValueType>
struct NormalFormData {
    storm::storage::BitVector const maybeStates;      // Those states that can be reached from initial without reaching a terminal state
    storm::storage::BitVector const terminalStates;   // Those states where we already know the probability to reach the condition and the target value
    storm::storage::BitVector const conditionStates;  // Those states where the condition holds almost surely (under all schedulers)
    storm::storage::BitVector const targetStates;     // Those states where the target holds almost surely (under all schedulers)
    storm::storage::BitVector const universalObservationFailureStates;    // Those states where the condition is not reachable (under all schedulers)
    storm::storage::BitVector const existentialObservationFailureStates;  // Those states s where a scheduler exists that (i) does not reach the condition from
                                                                          // s and (ii) acts optimal in all terminal states
    std::map<uint64_t, ValueType> const nonZeroTargetStateValues;         // The known non-zero target values. (default is zero)
    // There are three cases of terminal states:
    // 1. conditionStates: The condition holds, so the target value is the optimal probability to reach target from there
    // 2. targetStates: The target is reached, so the target value is the optimal probability to reach a condition from there.
    //                  The remaining probability mass is the probability of an observation failure
    // 3. states that can not reach the condition under any scheduler. The target value is zero.

    // TerminalStates is a superset of conditionStates and dom(nonZeroTargetStateValues).
    // For a terminalState that is not a conditionState, it is impossible to (reach the condition and not reach the target).

    std::vector<uint64_t> const
        schedulerChoicesForReachingTargetStates;  // Scheduler choices for reaching target states, used for constructing the resulting scheduler
    std::vector<uint64_t> const
        schedulerChoicesForReachingConditionStates;  // Scheduler choices for reaching condition states, used for constructing the resulting scheduler

    ValueType getTargetValue(uint64_t state) const {
        STORM_LOG_ASSERT(terminalStates.get(state), "Tried to get target value for non-terminal state");
        auto const it = nonZeroTargetStateValues.find(state);
        return it == nonZeroTargetStateValues.end() ? storm::utility::zero<ValueType>() : it->second;
    }

    ValueType failProbability(uint64_t state) const {
        STORM_LOG_ASSERT(terminalStates.get(state), "Tried to get fail probability for non-terminal state");
        STORM_LOG_ASSERT(!conditionStates.get(state), "Tried to get fail probability for a condition state");
        // condition states have fail probability zero
        return storm::utility::one<ValueType>() - getTargetValue(state);
    }
};

template<typename ValueType>
NormalFormData<ValueType> obtainNormalForm(Environment const& env, storm::solver::OptimizationDirection const dir,
                                           storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                           storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& relevantStates,
                                           storm::storage::BitVector const& targetStates, storm::storage::BitVector const& conditionStates) {
    storm::storage::BitVector const allStates(transitionMatrix.getRowGroupCount(), true);
    auto extendedConditionStates =
        storm::utility::graph::performProb1A(transitionMatrix, transitionMatrix.getRowGroupIndices(), backwardTransitions, allStates, conditionStates);
    auto universalObservationFailureStates = storm::utility::graph::performProb0A(backwardTransitions, allStates, extendedConditionStates);
    std::map<uint64_t, ValueType> nonZeroTargetStateValues;
    auto const extendedTargetStates =
        storm::utility::graph::performProb1A(transitionMatrix, transitionMatrix.getRowGroupIndices(), backwardTransitions, allStates, targetStates);
    auto const targetAndNotCondFailStates = extendedTargetStates & ~(extendedConditionStates | universalObservationFailureStates);

    std::vector<uint64_t> schedulerChoicesForReachingTargetStates;
    std::vector<uint64_t> schedulerChoicesForReachingConditionStates;

    auto schedulerForTargetStates = computeReachabilityProbabilitiesAndScheduler(env, nonZeroTargetStateValues, dir, transitionMatrix, extendedConditionStates,
                                                                                 allStates, extendedTargetStates);
    schedulerChoicesForReachingTargetStates = std::vector<uint64_t>(transitionMatrix.getRowGroupCount(), 0);
    if (schedulerForTargetStates) {
        for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
            schedulerChoicesForReachingTargetStates[state] = schedulerForTargetStates->getChoice(state).getDeterministicChoice();
        }
    }

    auto schedulerForConditionStates = computeReachabilityProbabilitiesAndScheduler(env, nonZeroTargetStateValues, dir, transitionMatrix,
                                                                                    targetAndNotCondFailStates, allStates, extendedConditionStates);
    schedulerChoicesForReachingConditionStates = std::vector<uint64_t>(transitionMatrix.getRowGroupCount(), 0);
    if (schedulerForConditionStates) {
        for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
            schedulerChoicesForReachingConditionStates[state] = schedulerForConditionStates->getChoice(state).getDeterministicChoice();
        }
    }

    // get states where the optimal policy reaches the condition with positive probability
    auto terminalStatesThatReachCondition = extendedConditionStates;
    for (auto state : targetAndNotCondFailStates) {
        if (nonZeroTargetStateValues.contains(state)) {
            terminalStatesThatReachCondition.set(state, true);
        }
    }

    // get the terminal states following the three cases described above
    auto terminalStates = extendedConditionStates | extendedTargetStates | universalObservationFailureStates;
    if (storm::solver::minimize(dir)) {
        // There can be target states from which (only) the *minimal* probability to reach a condition is zero.
        // For those states, the optimal policy is to enforce observation failure.
        // States that can only reach (target states with almost sure observation failure) or observation failure will be treated as terminal states with
        // targetValue zero and failProbability one.
        terminalStates |= storm::utility::graph::performProb0A(backwardTransitions, ~terminalStates, terminalStatesThatReachCondition);
    }

    auto nonTerminalStates = ~terminalStates;

    auto existentialObservationFailureStates = storm::utility::graph::performProb0E(transitionMatrix, transitionMatrix.getRowGroupIndices(),
                                                                                    backwardTransitions, nonTerminalStates, terminalStatesThatReachCondition);

    // Restrict non-terminal states to those that are still relevant
    nonTerminalStates &= storm::utility::graph::getReachableStates(transitionMatrix, relevantStates, nonTerminalStates, terminalStates);

    return NormalFormData<ValueType>{.maybeStates = std::move(nonTerminalStates),
                                     .terminalStates = std::move(terminalStates),
                                     .conditionStates = std::move(extendedConditionStates),
                                     .targetStates = std::move(extendedTargetStates),
                                     .universalObservationFailureStates = std::move(universalObservationFailureStates),
                                     .existentialObservationFailureStates = std::move(existentialObservationFailureStates),
                                     .nonZeroTargetStateValues = std::move(nonZeroTargetStateValues),
                                     .schedulerChoicesForReachingTargetStates = std::move(schedulerChoicesForReachingTargetStates),
                                     .schedulerChoicesForReachingConditionStates = std::move(schedulerChoicesForReachingConditionStates)};
}

// computes the scheduler that reaches the EC exits from the maybe states that were removed by EC elimination
template<typename ValueType, typename SolutionType = ValueType>
void finalizeSchedulerForMaybeStates(storm::storage::Scheduler<SolutionType>& scheduler, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                     storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& maybeStates,
                                     storm::storage::BitVector const& maybeStatesWithoutChoice, storm::storage::BitVector const& maybeStatesWithChoice,
                                     std::vector<uint64_t> const& stateToFinalEc, NormalFormData<ValueType> const& normalForm, uint64_t initialComponentIndex,
                                     storm::storage::BitVector const& initialComponentExitRows, uint64_t chosenInitialComponentExitState,
                                     uint64_t chosenInitialComponentExit) {
    // Compute the EC stay choices for the states in maybeStatesWithChoice
    storm::storage::BitVector ecStayChoices(transitionMatrix.getRowCount(), false);
    storm::storage::BitVector initialComponentStates(transitionMatrix.getRowGroupCount(), false);

    // compute initial component states and all choices that stay within a given EC
    for (auto state : maybeStates) {
        auto ecIndex = stateToFinalEc[state];
        if (ecIndex == initialComponentIndex) {
            initialComponentStates.set(state, true);
            continue;  // state part of the initial component
        } else if (ecIndex == std::numeric_limits<uint64_t>::max()) {
            continue;
        }
        for (auto choiceIndex : transitionMatrix.getRowGroupIndices(state)) {
            bool isEcStayChoice = true;
            for (auto const& entry : transitionMatrix.getRow(choiceIndex)) {
                auto targetState = entry.getColumn();
                if (stateToFinalEc[targetState] != ecIndex) {
                    isEcStayChoice = false;
                    break;
                }
            }
            if (isEcStayChoice) {
                ecStayChoices.set(choiceIndex, true);
            }
        }
    }

    // fill choices for ECs that reach the chosen EC exit
    auto const maybeNonICStatesWithoutChoice = maybeStatesWithoutChoice & ~initialComponentStates;
    storm::utility::graph::computeSchedulerProb1E(maybeNonICStatesWithoutChoice, transitionMatrix, backwardTransitions, maybeStates, maybeStatesWithChoice,
                                                  scheduler, ecStayChoices);

    // collect all choices from the initial component states and the choices that were selected by the scheduler so far
    auto const condOrTargetStates = normalForm.conditionStates | normalForm.targetStates;
    auto const rowGroups = transitionMatrix.getRowGroupIndices();
    storm::storage::BitVector allowedChoices(transitionMatrix.getRowCount(), false);
    auto const rowGroupCount = transitionMatrix.getRowGroupCount();
    for (uint64_t state = 0; state < rowGroupCount; ++state) {
        if (scheduler.isChoiceSelected(state)) {
            auto choiceIndex = scheduler.getChoice(state).getDeterministicChoice();
            allowedChoices.set(rowGroups[state] + choiceIndex, true);
        } else if (initialComponentStates.get(state) || condOrTargetStates.get(state)) {
            for (auto choiceIndex : transitionMatrix.getRowGroupIndices(state)) {
                allowedChoices.set(choiceIndex, true);
            }
        }
    }

    auto const transposedMatrixWithGroups = transitionMatrix.transpose(false);

    // dfs to find which choices in initial component states lead to condOrTargetStates
    storm::storage::BitVector choicesThatCanVisitCondOrTargetStates(transitionMatrix.getRowCount(), false);
    std::stack<uint64_t> toProcess;
    for (auto state : condOrTargetStates) {
        toProcess.push(state);
    }
    storm::storage::BitVector visitedStates(rowGroupCount, false);
    visitedStates = condOrTargetStates;
    while (!toProcess.empty()) {
        auto currentState = toProcess.top();
        toProcess.pop();
        for (auto const& entry : transposedMatrixWithGroups.getRow(currentState)) {
            auto predecessorChoiceIndex = entry.getColumn();
            if (!allowedChoices.get(predecessorChoiceIndex) || choicesThatCanVisitCondOrTargetStates.get(predecessorChoiceIndex)) {
                continue;
            }
            choicesThatCanVisitCondOrTargetStates.set(predecessorChoiceIndex, true);
            uint64_t predecessorState = 0;
            for (uint64_t state = 0; state < rowGroupCount; ++state) {
                if (predecessorChoiceIndex < rowGroups[state + 1]) {
                    predecessorState = state;
                    break;
                }
            }
            if (!visitedStates.get(predecessorState)) {
                visitedStates.set(predecessorState, true);
                toProcess.push(predecessorState);
            }
        }
    }

    // we want to disallow taking initial component exits that can lead to a condition or target state, beside the one exit that was chosen
    storm::storage::BitVector disallowedInitialComponentExits = initialComponentExitRows & choicesThatCanVisitCondOrTargetStates;
    disallowedInitialComponentExits.set(chosenInitialComponentExit, false);

    storm::storage::BitVector choicesAllowedForInitialComponent = allowedChoices & ~disallowedInitialComponentExits;
    storm::storage::BitVector exitStateBitvector(transitionMatrix.getRowGroupCount(), false);
    exitStateBitvector.set(chosenInitialComponentExitState, true);

    storm::utility::graph::computeSchedulerProbGreater0E(transitionMatrix, backwardTransitions, initialComponentStates, exitStateBitvector, scheduler,
                                                         choicesAllowedForInitialComponent);
}

template<typename ValueType, typename SolutionType = ValueType>
struct ResultReturnType {
    ResultReturnType(ValueType initialStateValue, std::unique_ptr<storm::storage::Scheduler<ValueType>>&& scheduler = nullptr)
        : initialStateValue(initialStateValue), scheduler(std::move(scheduler)) {
        // Intentionally left empty.
    }

    bool hasScheduler() const {
        return static_cast<bool>(scheduler);
    }

    ValueType initialStateValue;
    std::unique_ptr<storm::storage::Scheduler<SolutionType>> scheduler;
};

/*!
 * Uses the restart method by Baier et al.
// @see doi.org/10.1007/978-3-642-54862-8_43
 */
template<typename ValueType, typename SolutionType = ValueType>
typename internal::ResultReturnType<ValueType> computeViaRestartMethod(Environment const& env, uint64_t const initialState,
                                                                       storm::solver::OptimizationDirection const dir,
                                                                       storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                       storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                                       NormalFormData<ValueType> const& normalForm) {
    auto const& maybeStates = normalForm.maybeStates;
    auto const stateToMatrixIndexMap = maybeStates.getNumberOfSetBitsBeforeIndices();
    auto const numMaybeStates = maybeStates.getNumberOfSetBits();
    auto const numMaybeChoices = transitionMatrix.getNumRowsInRowGroups(maybeStates);

    // Build the transitions that include a backwards loop to the initial state
    storm::storage::SparseMatrixBuilder<ValueType> matrixBuilder(numMaybeChoices, numMaybeStates, 0, true, true, numMaybeStates);
    std::vector<ValueType> rowValues;
    storm::storage::BitVector rowsWithSum1(numMaybeChoices, true);
    rowValues.reserve(numMaybeChoices);
    uint64_t currentRow = 0;
    for (auto state : maybeStates) {
        matrixBuilder.newRowGroup(currentRow);
        for (auto origRowIndex : transitionMatrix.getRowGroupIndices(state)) {
            // We make two passes over the successors. First, we find out the reset probabilities and target probabilities
            // Then, we insert the matrix entries in the correct order
            // This two-phase approach is to avoid a costly out-of-order insertion into the matrix
            ValueType targetProbability = storm::utility::zero<ValueType>();
            ValueType restartProbability = storm::utility::zero<ValueType>();
            bool rowSumIsLess1 = false;
            for (auto const& entry : transitionMatrix.getRow(origRowIndex)) {
                if (normalForm.terminalStates.get(entry.getColumn())) {
                    ValueType const targetValue = normalForm.getTargetValue(entry.getColumn());
                    targetProbability += targetValue * entry.getValue();
                    if (normalForm.conditionStates.get(entry.getColumn())) {
                        rowSumIsLess1 = true;
                    } else {
                        if (!storm::utility::isZero(targetValue)) {
                            rowSumIsLess1 = true;
                        }
                        restartProbability += entry.getValue() * normalForm.failProbability(entry.getColumn());
                    }
                }
            }
            if (rowSumIsLess1) {
                rowsWithSum1.set(currentRow, false);
            }
            rowValues.push_back(targetProbability);
            bool addRestartTransition = !storm::utility::isZero(restartProbability);
            for (auto const& entry : transitionMatrix.getRow(origRowIndex)) {
                // Insert backloop probability if we haven't done so yet and are past the initial state index
                // This is to avoid a costly out-of-order insertion into the matrix
                if (addRestartTransition && entry.getColumn() > initialState) {
                    matrixBuilder.addNextValue(currentRow, stateToMatrixIndexMap[initialState], restartProbability);
                    addRestartTransition = false;
                }
                if (maybeStates.get(entry.getColumn())) {
                    matrixBuilder.addNextValue(currentRow, stateToMatrixIndexMap[entry.getColumn()], entry.getValue());
                }
            }
            // Add the backloop if we haven't done this already
            if (addRestartTransition) {
                matrixBuilder.addNextValue(currentRow, stateToMatrixIndexMap[initialState], restartProbability);
            }
            ++currentRow;
        }
    }

    auto maybeMatrix = matrixBuilder.build();
    auto matrix = storm::storage::SparseMatrix<ValueType>(maybeMatrix);
    auto initStateInMatrix = stateToMatrixIndexMap[initialState];

    // Eliminate end components in two phases
    // First, we catch all end components that do not contain the initial state. It is possible to stay in those ECs forever
    // without reaching the condition. This is reflected by a backloop to the initial state.
    storm::storage::BitVector selectedStatesInMatrix(numMaybeStates, true);
    selectedStatesInMatrix.set(initStateInMatrix, false);
    auto ecElimResult1 = eliminateEndComponents(selectedStatesInMatrix, true, initStateInMatrix, matrix, initStateInMatrix, rowsWithSum1, rowValues);
    // Second, eliminate the remaining ECs. These must involve the initial state and might have been introduced in the previous step.
    // A policy selecting such an EC must reach the condition with probability zero and is thus invalid.
    selectedStatesInMatrix.set(initStateInMatrix, true);
    auto ecElimResult2 = eliminateEndComponents(selectedStatesInMatrix, false, std::nullopt, matrix, initStateInMatrix, rowsWithSum1, rowValues);

    STORM_LOG_INFO("Processed model has " << matrix.getRowGroupCount() << " states and " << matrix.getRowGroupCount() << " choices and "
                                          << matrix.getEntryCount() << " transitions.");
    // Finally, solve the equation system
    auto result = solveMinMaxEquationSystem(env, matrix, rowValues, rowsWithSum1, dir, initStateInMatrix);

    storm::storage::BitVector initialComponentExitRows(transitionMatrix.getRowCount(), false);
    for (auto rowIndex : matrix.getRowGroupIndices(initStateInMatrix)) {
        uint64_t originalRowIndex = rowIndex;
        if (ecElimResult2.has_value()) {
            originalRowIndex = ecElimResult2->newToOldRowMapping[originalRowIndex];
        }
        if (ecElimResult1.has_value()) {
            originalRowIndex = ecElimResult1->newToOldRowMapping[originalRowIndex];
        }

        uint64_t originalState;
        uint64_t originalChoice;
        auto const rowGroups = maybeMatrix.getRowGroupIndices();
        for (originalState = 0; originalState < maybeMatrix.getRowGroupCount(); ++originalState) {
            auto const firstRowStateIndex = rowGroups[originalState + 1];
            if (firstRowStateIndex > originalRowIndex) {
                originalChoice = originalRowIndex - rowGroups[originalState];
                break;
            }
        }

        uint64_t index = maybeStates.getNextSetIndex(0);
        for (uint64_t s = 0; s < originalState; ++s) {
            index = maybeStates.getNextSetIndex(index + 1);
        }

        originalState = index;

        initialComponentExitRows.set(transitionMatrix.getRowGroupIndices()[originalState] + originalChoice, true);
    }

    // std::vector<uint64_t> finalSchedulerChoices(transitionMatrix.getRowGroupCount(), -1);
    storm::storage::BitVector maybeStatesWithChoice(maybeStates.size(), false);
    std::unique_ptr<storm::storage::Scheduler<SolutionType>> scheduler;
    uint64_t chosenInitialComponentExitState;
    uint64_t chosenInitialComponentExit;
    scheduler = std::make_unique<storm::storage::Scheduler<SolutionType>>(transitionMatrix.getRowGroupCount());
    // std::vector<uint64_t> ecsToExits(matrix.getRowGroupCount(), 0);
    if (result.hasScheduler()) {
        uint64_t state = 0;
        for (auto& choice : *result.scheduler) {
            uint64_t firstRowIndex = matrix.getRowGroupIndices()[state];
            uint64_t originalChoice = firstRowIndex + choice;
            if (ecElimResult2.has_value()) {
                originalChoice = ecElimResult2->newToOldRowMapping[originalChoice];
            }
            if (ecElimResult1.has_value()) {
                originalChoice = ecElimResult1->newToOldRowMapping[originalChoice];
            }

            uint64_t originalState;
            auto const rowGroups = maybeMatrix.getRowGroupIndices();
            for (originalState = 0; originalState < maybeMatrix.getRowGroupCount(); ++originalState) {
                auto const firstRowStateIndex = rowGroups[originalState + 1];
                if (firstRowStateIndex > originalChoice) {
                    originalChoice = originalChoice - rowGroups[originalState];
                    break;
                }
            }

            uint64_t index = maybeStates.getNextSetIndex(0);
            for (uint64_t s = 0; s < originalState; ++s) {
                index = maybeStates.getNextSetIndex(index + 1);
            }

            originalState = index;
            // finalSchedulerChoices[originalState] = originalChoice;
            scheduler->setChoice(originalChoice, originalState);
            maybeStatesWithChoice.set(originalState, true);
            // ecsToExits[state] = originalState;
            if (state == initStateInMatrix) {
                chosenInitialComponentExitState = originalState;
                chosenInitialComponentExit = transitionMatrix.getRowGroupIndices()[originalState] + originalChoice;
            }
            ++state;
        }
    }

    std::vector<uint64_t> stateToFinalEc(transitionMatrix.getRowGroupCount(), std::numeric_limits<uint64_t>::max());
    uint64_t state = 0;
    for (auto s : maybeStates) {
        auto mappedState = state;
        mappedState = ecElimResult1.has_value() ? ecElimResult1->oldToNewStateMapping[mappedState] : mappedState;
        mappedState = ecElimResult2.has_value() ? ecElimResult2->oldToNewStateMapping[mappedState] : mappedState;
        stateToFinalEc[s] = mappedState;
        state++;
    }

    auto const maybeStatesWithoutChoice = maybeStates & ~maybeStatesWithChoice;
    finalizeSchedulerForMaybeStates(*scheduler, transitionMatrix, backwardTransitions, maybeStates, maybeStatesWithoutChoice, maybeStatesWithChoice,
                                    stateToFinalEc, normalForm, initStateInMatrix, initialComponentExitRows, chosenInitialComponentExitState,
                                    chosenInitialComponentExit);

    auto finalResult = ResultReturnType<ValueType>(result.initialStateValue, std::move(scheduler));

    return finalResult;
}

/*!
 * A helper class that computes (weighted) reachability probabilities for a given MDP in normal form.
 * @tparam ValueType
 * @tparam SolutionType
 */
template<typename ValueType, typename SolutionType = ValueType>
class WeightedReachabilityHelper {
   public:
    WeightedReachabilityHelper(uint64_t const initialState, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                               NormalFormData<ValueType> const& normalForm) {
        // Determine rowgroups (states) and rows (choices) of the submatrix
        auto subMatrixRowGroups = normalForm.maybeStates;
        // Identify and eliminate the initial component to enforce that it is eventually exited
        // The initial component is the largest subset of maybestates C such that
        // (i) the initial state is contained in C
        // (ii) each state in C can be reached from the initial state while only playing actions that stay inside C or observation failure and
        // (iii) for each state in C except the initial state there is a policy that almost surely reaches an observation failure
        // An optimal scheduler can intuitively pick the best exiting action of C and enforce that all paths that satisfy the condition exit C through that
        // action. By eliminating the initial component, we ensure that only policies that actually exit C are considered. The remaining policies have
        // probability zero of satisfying the condition.
        initialComponentExitRows = storm::storage::BitVector(transitionMatrix.getRowCount(), false);
        subMatrixRowGroups.set(initialState, false);  // temporarily unset initial state
        std::vector<uint64_t> dfsStack = {initialState};
        while (!dfsStack.empty()) {
            auto const state = dfsStack.back();
            dfsStack.pop_back();
            for (auto rowIndex : transitionMatrix.getRowGroupIndices(state)) {
                auto const row = transitionMatrix.getRow(rowIndex);
                if (std::all_of(row.begin(), row.end(),
                                [&normalForm](auto const& entry) { return normalForm.existentialObservationFailureStates.get(entry.getColumn()); })) {
                    for (auto const& entry : row) {
                        auto const& successor = entry.getColumn();
                        if (subMatrixRowGroups.get(successor)) {
                            subMatrixRowGroups.set(successor, false);
                            dfsStack.push_back(successor);
                        }
                    }
                } else {
                    initialComponentExitRows.set(rowIndex, true);
                }
            }
        }
        auto const numSubmatrixRows = transitionMatrix.getNumRowsInRowGroups(subMatrixRowGroups) + initialComponentExitRows.getNumberOfSetBits();
        subMatrixRowGroups.set(initialState, true);  // set initial state again, as single representative state for the initial component
        auto const numSubmatrixRowGroups = subMatrixRowGroups.getNumberOfSetBits();

        // state index mapping and initial state
        auto stateToMatrixIndexMap = subMatrixRowGroups.getNumberOfSetBitsBeforeIndices();
        initialStateInSubmatrix = stateToMatrixIndexMap[initialState];
        auto const eliminatedInitialComponentStates = normalForm.maybeStates & ~subMatrixRowGroups;
        for (auto state : eliminatedInitialComponentStates) {
            stateToMatrixIndexMap[state] = initialStateInSubmatrix;  // map all eliminated states to the initial state
        }

        // build matrix, rows that sum up to 1, target values, condition values
        storm::storage::SparseMatrixBuilder<ValueType> matrixBuilder(numSubmatrixRows, numSubmatrixRowGroups, 0, true, true, numSubmatrixRowGroups);
        rowsWithSum1 = storm::storage::BitVector(numSubmatrixRows, true);
        targetRowValues.reserve(numSubmatrixRows);
        conditionRowValues.reserve(numSubmatrixRows);
        uint64_t currentRow = 0;
        for (auto state : subMatrixRowGroups) {
            matrixBuilder.newRowGroup(currentRow);

            // Put the row processing into a lambda for avoiding code duplications
            auto processRow = [&](uint64_t origRowIndex) {
                // We make two passes. First, we find out the probability to reach an eliminated initial component state
                ValueType const eliminatedInitialComponentProbability = transitionMatrix.getConstrainedRowSum(origRowIndex, eliminatedInitialComponentStates);
                // Second, we insert the submatrix entries and find out the target and condition probabilities for this row
                ValueType targetProbability = storm::utility::zero<ValueType>();
                ValueType conditionProbability = storm::utility::zero<ValueType>();
                bool rowSumIsLess1 = false;
                bool initialStateEntryInserted = false;
                for (auto const& entry : transitionMatrix.getRow(origRowIndex)) {
                    if (normalForm.terminalStates.get(entry.getColumn())) {
                        STORM_LOG_ASSERT(!storm::utility::isZero(entry.getValue()), "Transition probability must be non-zero");
                        rowSumIsLess1 = true;
                        ValueType const scaledTargetValue = normalForm.getTargetValue(entry.getColumn()) * entry.getValue();
                        targetProbability += scaledTargetValue;
                        if (normalForm.conditionStates.get(entry.getColumn())) {
                            conditionProbability += entry.getValue();  // conditionValue of successor is 1
                        } else {
                            conditionProbability += scaledTargetValue;  // for terminal, non-condition states, the condition value equals the target value
                        }
                    } else if (!eliminatedInitialComponentStates.get(entry.getColumn())) {
                        auto const columnIndex = stateToMatrixIndexMap[entry.getColumn()];
                        if (!initialStateEntryInserted && columnIndex >= initialStateInSubmatrix) {
                            if (columnIndex == initialStateInSubmatrix) {
                                matrixBuilder.addNextValue(currentRow, initialStateInSubmatrix, eliminatedInitialComponentProbability + entry.getValue());
                            } else {
                                matrixBuilder.addNextValue(currentRow, initialStateInSubmatrix, eliminatedInitialComponentProbability);
                                matrixBuilder.addNextValue(currentRow, columnIndex, entry.getValue());
                            }
                            initialStateEntryInserted = true;
                        } else {
                            matrixBuilder.addNextValue(currentRow, columnIndex, entry.getValue());
                        }
                    }
                }
                if (rowSumIsLess1) {
                    rowsWithSum1.set(currentRow, false);
                }
                targetRowValues.push_back(targetProbability);
                conditionRowValues.push_back(conditionProbability);
                ++currentRow;
            };
            // invoke the lambda
            if (state == initialState) {
                for (auto origRowIndex : initialComponentExitRows) {
                    processRow(origRowIndex);
                    initialComponentExitToOriginalRow.push_back(origRowIndex);
                }
            } else {
                for (auto origRowIndex : transitionMatrix.getRowGroupIndices(state)) {
                    processRow(origRowIndex);
                }
            }
        }
        fullSubmatrix = matrixBuilder.build();
        submatrix = storm::storage::SparseMatrix<ValueType>(fullSubmatrix);

        //  eliminate ECs if present. We already checked that the initial state can not yield observation failure, so it cannot be part of an EC.
        //  For all remaining ECs, staying in an EC forever is reflected by collecting a value of zero for both, target and condition
        storm::storage::BitVector allExceptInit(numSubmatrixRowGroups, true);
        allExceptInit.set(initialStateInSubmatrix, false);
        ecResult = eliminateEndComponents<ValueType>(allExceptInit, true, std::nullopt, submatrix, initialStateInSubmatrix, rowsWithSum1, targetRowValues,
                                                     conditionRowValues);
        STORM_LOG_INFO("Processed model has " << submatrix.getRowGroupCount() << " states and " << submatrix.getRowGroupCount() << " choices and "
                                              << submatrix.getEntryCount() << " transitions.");

        // if (ecElimResult.has_value()) {
        //     ecResult = std::make_unique<typename
        //     storm::transformer::EndComponentEliminator<ValueType>::EndComponentEliminatorReturnType>(std::move(*ecElimResult));
        // }

        stateToFinalEc.resize(transitionMatrix.getRowGroupCount(), -1);
        auto maybeStatesNotInSubmatrix = normalForm.maybeStates & ~subMatrixRowGroups;

        for (auto state : maybeStatesNotInSubmatrix) {
            stateToFinalEc[state] = 0;
        }

        uint64_t state = 0;
        for (auto s : subMatrixRowGroups) {
            uint64_t mappedState = state;
            mappedState = ecResult.has_value() ? ecResult->oldToNewStateMapping[mappedState] : mappedState;
            stateToFinalEc[s] = mappedState;
            state++;
        }
    }

    internal::SolverResult<ValueType> computeWeightedDiff(storm::Environment const& env, storm::OptimizationDirection const dir, ValueType const& targetWeight,
                                                          ValueType const& conditionWeight) const {
        auto rowValues = createScaledVector(targetWeight, targetRowValues, conditionWeight, conditionRowValues);

        // Initialize the solution vector.
        std::vector<SolutionType> x(submatrix.getRowGroupCount(), storm::utility::zero<ValueType>());

        // Set up the solver.
        auto solver = storm::solver::GeneralMinMaxLinearEquationSolverFactory<ValueType, SolutionType>().create(env, submatrix);
        solver->setOptimizationDirection(dir);
        solver->setRequirementsChecked();
        solver->setHasUniqueSolution(true);
        solver->setHasNoEndComponents(true);
        solver->setLowerBound(-storm::utility::one<ValueType>());
        solver->setUpperBound(storm::utility::one<ValueType>());
        solver->setTrackScheduler(true);

        // Solve the corresponding system of equations.
        solver->solveEquations(env, x, rowValues);

        SolverResult<SolutionType> result(x[initialStateInSubmatrix]);
        result.scheduler = std::move(solver->getSchedulerChoices());

        return result;
    }

    auto getInternalInitialState() const {
        return initialStateInSubmatrix;
    }

    void evaluateScheduler(storm::Environment const& env, std::vector<uint64_t>& scheduler, std::vector<SolutionType>& targetResults,
                           std::vector<SolutionType>& conditionResults) const {
        if (scheduler.empty()) {
            scheduler.resize(submatrix.getRowGroupCount(), 0);
        }
        if (targetResults.empty()) {
            targetResults.resize(submatrix.getRowGroupCount(), storm::utility::zero<ValueType>());
        }
        if (conditionResults.empty()) {
            conditionResults.resize(submatrix.getRowGroupCount(), storm::utility::zero<ValueType>());
        }
        // apply the scheduler
        storm::solver::GeneralLinearEquationSolverFactory<ValueType> factory;
        bool const convertToEquationSystem = factory.getEquationProblemFormat(env) == storm::solver::LinearEquationSolverProblemFormat::EquationSystem;
        auto scheduledMatrix = submatrix.selectRowsFromRowGroups(scheduler, convertToEquationSystem);
        if (convertToEquationSystem) {
            scheduledMatrix.convertToEquationSystem();
        }
        auto solver = factory.create(env, std::move(scheduledMatrix));
        solver->setBounds(storm::utility::zero<ValueType>(), storm::utility::one<ValueType>());
        solver->setCachingEnabled(true);

        std::vector<ValueType> subB(submatrix.getRowGroupCount());
        storm::utility::vector::selectVectorValues<ValueType>(subB, scheduler, submatrix.getRowGroupIndices(), targetRowValues);
        solver->solveEquations(env, targetResults, subB);

        storm::utility::vector::selectVectorValues<ValueType>(subB, scheduler, submatrix.getRowGroupIndices(), conditionRowValues);
        solver->solveEquations(env, conditionResults, subB);
    }

    template<OptimizationDirection Dir>
    bool improveScheduler(std::vector<uint64_t>& scheduler, ValueType const& lambda, std::vector<SolutionType> const& targetResults,
                          std::vector<SolutionType> const& conditionResults) {
        bool improved = false;
        for (uint64_t rowGroupIndex = 0; rowGroupIndex < scheduler.size(); ++rowGroupIndex) {
            storm::utility::Extremum<Dir, ValueType> groupValue;
            uint64_t optimalRowIndex{0};
            ValueType scheduledValue;
            for (auto rowIndex : submatrix.getRowGroupIndices(rowGroupIndex)) {
                ValueType rowValue = targetRowValues[rowIndex] - lambda * conditionRowValues[rowIndex];
                for (auto const& entry : submatrix.getRow(rowIndex)) {
                    rowValue += entry.getValue() * (targetResults[entry.getColumn()] - lambda * conditionResults[entry.getColumn()]);
                }
                if (rowIndex == scheduler[rowGroupIndex] + submatrix.getRowGroupIndices()[rowGroupIndex]) {
                    scheduledValue = rowValue;
                }
                if (groupValue &= rowValue) {
                    optimalRowIndex = rowIndex;
                }
            }
            if (scheduledValue != *groupValue) {
                scheduler[rowGroupIndex] = optimalRowIndex - submatrix.getRowGroupIndices()[rowGroupIndex];
                improved = true;
            }
        }
        return improved;
    }

    storm::storage::SparseMatrix<ValueType> submatrix;
    storm::storage::SparseMatrix<ValueType> fullSubmatrix;
    std::vector<uint64_t> stateToFinalEc;
    boost::optional<typename storm::transformer::EndComponentEliminator<ValueType>::EndComponentEliminatorReturnType> ecResult;
    std::vector<uint64_t> initialComponentExitToOriginalRow;
    storm::storage::BitVector initialComponentExitRows;

   private:
    std::vector<ValueType> createScaledVector(ValueType const& w1, std::vector<ValueType> const& v1, ValueType const& w2,
                                              std::vector<ValueType> const& v2) const {
        STORM_LOG_ASSERT(v1.size() == v2.size(), "Vector sizes must match");
        std::vector<ValueType> result;
        result.reserve(v1.size());
        for (size_t i = 0; i < v1.size(); ++i) {
            result.push_back(w1 * v1[i] + w2 * v2[i]);
        }
        return result;
    }

    storm::storage::BitVector rowsWithSum1;
    std::vector<ValueType> targetRowValues;
    std::vector<ValueType> conditionRowValues;
    uint64_t initialStateInSubmatrix;
};

enum class BisectionMethodBounds { Simple, Advanced };
template<typename ValueType, typename SolutionType = ValueType>
typename internal::ResultReturnType<ValueType> computeViaBisection(Environment const& env, BisectionMethodBounds boundOption, uint64_t const initialState,
                                                                   storm::solver::SolveGoal<ValueType, SolutionType> goal,
                                                                   storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                   storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                                   NormalFormData<ValueType> const& normalForm) {
    // We currently handle sound model checking incorrectly: we would need the actual lower/upper bounds of the weightedReachabilityHelper
    STORM_LOG_WARN_COND(!env.solver().isForceSoundness(),
                        "Bisection method does not adequately handle propagation of errors. Result is not necessarily sound.");
    SolutionType const precision = [&env, boundOption]() {
        if (storm::NumberTraits<SolutionType>::IsExact || env.solver().isForceExact()) {
            STORM_LOG_WARN_COND(storm::NumberTraits<SolutionType>::IsExact && boundOption == BisectionMethodBounds::Advanced,
                                "Selected bisection method with exact precision in a setting that might not terminate.");
            return storm::utility::zero<SolutionType>();
        } else {
            return storm::utility::convertNumber<SolutionType>(env.solver().minMax().getPrecision());
        }
    }();
    bool const relative = env.solver().minMax().getRelativeTerminationCriterion();

    WeightedReachabilityHelper wrh(initialState, transitionMatrix, normalForm);
    SolutionType pMin{storm::utility::zero<SolutionType>()};
    SolutionType pMax{storm::utility::one<SolutionType>()};

    if (boundOption == BisectionMethodBounds::Advanced) {
        auto pMinRes =
            wrh.computeWeightedDiff(env, storm::OptimizationDirection::Minimize, storm::utility::zero<ValueType>(), storm::utility::one<ValueType>());
        auto pMaxRes =
            wrh.computeWeightedDiff(env, storm::OptimizationDirection::Maximize, storm::utility::zero<ValueType>(), storm::utility::one<ValueType>());
        pMin = pMinRes.initialStateValue;
        pMax = pMaxRes.initialStateValue;
        STORM_LOG_TRACE("Conditioning event bounds:\n\t Lower bound: " << storm::utility::convertNumber<double>(pMin)
                                                                       << ",\n\t Upper bound: " << storm::utility::convertNumber<double>(pMax));
    }
    storm::utility::Extremum<storm::OptimizationDirection::Maximize, SolutionType> lowerBound = storm::utility::zero<ValueType>();
    storm::utility::Extremum<storm::OptimizationDirection::Minimize, SolutionType> upperBound = storm::utility::one<ValueType>();
    SolutionType middle;
    if (goal.isBounded()) {
        middle = goal.thresholdValue();
    } else {
        middle = (*lowerBound + *upperBound) / 2;
    }
    SolverResult<ValueType> res(0);
    for (uint64_t iterationCount = 1; true; ++iterationCount) {
        // evaluate the current middle
        res = wrh.computeWeightedDiff(env, goal.direction(), storm::utility::one<ValueType>(), -middle);
        SolutionType const middleValue = res.initialStateValue;
        // update the bounds and new middle value according to the bisection method
        if (boundOption == BisectionMethodBounds::Simple) {
            if (middleValue >= storm::utility::zero<ValueType>()) {
                lowerBound &= middle;
            }
            if (middleValue <= storm::utility::zero<ValueType>()) {
                upperBound &= middle;
            }
            middle = (*lowerBound + *upperBound) / 2;  // update middle to the average of the bounds
        } else {
            STORM_LOG_ASSERT(boundOption == BisectionMethodBounds::Advanced, "Unknown bisection method bounds");
            if (middleValue >= storm::utility::zero<ValueType>()) {
                lowerBound &= middle + (middleValue / pMax);
                upperBound &= middle + (middleValue / pMin);
            }
            if (middleValue <= storm::utility::zero<ValueType>()) {
                lowerBound &= middle + (middleValue / pMin);
                upperBound &= middle + (middleValue / pMax);
            }
            // update middle to the average of the bounds, but scale it according to the middle value (which is in [-1,1])
            middle = *lowerBound + (storm::utility::one<SolutionType>() + middleValue) * (*upperBound - *lowerBound) / 2;

            if (!storm::NumberTraits<SolutionType>::IsExact && storm::utility::isAlmostZero(*upperBound - *lowerBound)) {
                if (*lowerBound > *upperBound) {
                    std::swap(*lowerBound, *upperBound);
                }
                STORM_LOG_WARN("Precision of non-exact type reached during bisection method. Result might be inaccurate.");
            } else {
                STORM_LOG_ASSERT(middle >= *lowerBound && middle <= *upperBound, "Bisection method bounds are inconsistent.");
            }
        }
        // check for convergence
        SolutionType const boundDiff = *upperBound - *lowerBound;
        STORM_LOG_TRACE("Iteration #" << iterationCount << ":\n\t Lower bound:      " << *lowerBound << ",\n\t Upper bound:      " << *upperBound
                                      << ",\n\t Difference:       " << boundDiff << ",\n\t Middle val:       " << middleValue
                                      << ",\n\t Difference bound: " << (relative ? (precision * *lowerBound) : precision) << ".");
        if (goal.isBounded()) {
            STORM_LOG_TRACE("Using threshold " << storm::utility::convertNumber<double>(goal.thresholdValue()) << " with comparison "
                                               << (goal.boundIsALowerBound() ? (goal.boundIsStrict() ? ">" : ">=") : (goal.boundIsStrict() ? "<" : "<="))
                                               << ".");
        }
        if (boundDiff <= (relative ? (precision * *lowerBound) : precision)) {
            STORM_LOG_INFO("Bisection method converged after " << iterationCount << " iterations. Difference is "
                                                               << std::setprecision(std::numeric_limits<double>::digits10)
                                                               << storm::utility::convertNumber<double>(boundDiff) << ".");
            break;
        }
        // Check if bounds are fully below or above threshold
        if (goal.isBounded() && (*upperBound <= goal.thresholdValue() || (*lowerBound >= goal.thresholdValue()))) {
            STORM_LOG_INFO("Bisection method determined result after " << iterationCount << " iterations. Found bounds are ["
                                                                       << storm::utility::convertNumber<double>(*lowerBound) << ", "
                                                                       << storm::utility::convertNumber<double>(*upperBound) << "], threshold is "
                                                                       << storm::utility::convertNumber<double>(goal.thresholdValue()) << ".");
            break;
        }
        // check for early termination
        if (storm::utility::resources::isTerminate()) {
            STORM_LOG_WARN("Bisection solver aborted after " << iterationCount << "iterations. Bound difference is "
                                                             << storm::utility::convertNumber<double>(boundDiff) << ".");
            break;
        }
        // process the middle value for the next iteration
        // This sets the middle value to a rational number with the smallest enumerator/denominator that is still within the bounds
        // With close bounds this can lead to the middle being set to exactly the lower or upper bound, thus allowing for an exact answer.
        if constexpr (storm::NumberTraits<SolutionType>::IsExact) {
            // find a rational number with a concise representation close to middle and within the bounds
            auto const exactMiddle = middle;

            // Find number of digits - 1. Method using log10 does not work since that uses doubles internally.
            auto numDigits = storm::utility::numDigits<SolutionType>(*upperBound - *lowerBound) - 1;

            do {
                ++numDigits;
                middle = storm::utility::kwek_mehlhorn::sharpen<SolutionType, SolutionType>(numDigits, exactMiddle);
            } while (middle <= *lowerBound || middle >= *upperBound);
        }
        // Since above code never sets 'middle' to exactly zero or one, we check if that could be necessary after a couple of iterations
        if (iterationCount == 8) {  // 8 is just a heuristic value, it could be any number
            if (storm::utility::isZero(*lowerBound)) {
                middle = storm::utility::zero<SolutionType>();
            } else if (storm::utility::isOne(*upperBound)) {
                middle = storm::utility::one<SolutionType>();
            }
        }
    }

    storm::storage::BitVector maybeStatesWithChoice(normalForm.maybeStates.size(), false);
    std::unique_ptr<storm::storage::Scheduler<SolutionType>> scheduler;
    uint64_t chosenInitialComponentExitState;
    uint64_t chosenInitialComponentExit;
    scheduler = std::make_unique<storm::storage::Scheduler<SolutionType>>(transitionMatrix.getRowGroupCount());
    if (res.hasScheduler()) {
        uint64_t state = 0;
        for (auto& choice : *res.scheduler) {
            uint64_t originalChoice;
            uint64_t originalState;

            if (state == wrh.getInternalInitialState()) {
                originalChoice = wrh.initialComponentExitToOriginalRow[choice];

                auto const rowGroups = transitionMatrix.getRowGroupIndices();
                for (originalState = 0; originalState < transitionMatrix.getRowGroupCount(); ++originalState) {
                    auto const firstRowStateIndex = rowGroups[originalState + 1];
                    if (firstRowStateIndex > originalChoice) {
                        originalChoice = originalChoice - rowGroups[originalState];
                        break;
                    }
                }

                scheduler->setChoice(originalChoice, originalState);
                maybeStatesWithChoice.set(originalState, true);
                chosenInitialComponentExitState = originalState;
                chosenInitialComponentExit = transitionMatrix.getRowGroupIndices()[originalState] + originalChoice;
                ++state;
                continue;
            }

            uint64_t firstRowIndex = wrh.submatrix.getRowGroupIndices()[state];
            originalChoice = firstRowIndex + choice;
            if (wrh.ecResult.has_value()) {
                originalChoice = wrh.ecResult->newToOldRowMapping[originalChoice];
            }

            auto const rowGroups = wrh.fullSubmatrix.getRowGroupIndices();
            for (originalState = 0; originalState < wrh.fullSubmatrix.getRowGroupCount(); ++originalState) {
                auto const firstRowStateIndex = rowGroups[originalState + 1];
                if (firstRowStateIndex > originalChoice) {
                    originalChoice = originalChoice - rowGroups[originalState];
                    break;
                }
            }

            uint64_t index = normalForm.maybeStates.getNextSetIndex(0);
            for (uint64_t s = 0; s < originalState; ++s) {
                index = normalForm.maybeStates.getNextSetIndex(index + 1);
            }

            originalState = index;
            scheduler->setChoice(originalChoice, originalState);
            maybeStatesWithChoice.set(originalState, true);
            ++state;
        }
    }

    auto const maybeStatesWithoutChoice = normalForm.maybeStates & ~maybeStatesWithChoice;
    finalizeSchedulerForMaybeStates(*scheduler, transitionMatrix, backwardTransitions, normalForm.maybeStates, maybeStatesWithoutChoice, maybeStatesWithChoice,
                                    wrh.stateToFinalEc, normalForm, wrh.getInternalInitialState(), wrh.initialComponentExitRows,
                                    chosenInitialComponentExitState, chosenInitialComponentExit);

    auto finalResult = ResultReturnType<ValueType>((*lowerBound + *upperBound) / 2, std::move(scheduler));

    return finalResult;
}

template<typename ValueType, typename SolutionType = ValueType>
SolutionType computeViaPolicyIteration(Environment const& env, uint64_t const initialState, storm::solver::OptimizationDirection const dir,
                                       storm::storage::SparseMatrix<ValueType> const& transitionMatrix, NormalFormData<ValueType> const& normalForm) {
    WeightedReachabilityHelper wrh(initialState, transitionMatrix, normalForm);

    std::vector<uint64_t> scheduler;
    std::vector<SolutionType> targetResults, conditionResults;
    for (uint64_t iterationCount = 1; true; ++iterationCount) {
        wrh.evaluateScheduler(env, scheduler, targetResults, conditionResults);
        STORM_LOG_WARN_COND(
            targetResults[wrh.getInternalInitialState()] <= conditionResults[wrh.getInternalInitialState()],
            "Potential numerical issues: the probability to reach the target is greater than the probability to reach the condition. Difference is "
                << (storm::utility::convertNumber<double, ValueType>(targetResults[wrh.getInternalInitialState()] -
                                                                     conditionResults[wrh.getInternalInitialState()]))
                << ".");
        ValueType const lambda = storm::utility::isZero(conditionResults[wrh.getInternalInitialState()])
                                     ? storm::utility::zero<ValueType>()
                                     : ValueType(targetResults[wrh.getInternalInitialState()] / conditionResults[wrh.getInternalInitialState()]);
        bool schedulerChanged{false};
        if (storm::solver::minimize(dir)) {
            schedulerChanged = wrh.template improveScheduler<storm::OptimizationDirection::Minimize>(scheduler, lambda, targetResults, conditionResults);
        } else {
            schedulerChanged = wrh.template improveScheduler<storm::OptimizationDirection::Maximize>(scheduler, lambda, targetResults, conditionResults);
        }
        if (!schedulerChanged) {
            STORM_LOG_INFO("Policy iteration for conditional probabilities converged after " << iterationCount << " iterations.");
            return lambda;
        }
        if (storm::utility::resources::isTerminate()) {
            STORM_LOG_WARN("Policy iteration for conditional probabilities converged aborted after " << iterationCount << "iterations.");
            return lambda;
        }
    }
}

template<typename ValueType, typename SolutionType = ValueType>
std::optional<SolutionType> handleTrivialCases(uint64_t const initialState, NormalFormData<ValueType> const& normalForm) {
    if (normalForm.terminalStates.get(initialState)) {
        STORM_LOG_DEBUG("Initial state is terminal.");
        if (normalForm.conditionStates.get(initialState)) {
            return normalForm.getTargetValue(initialState);  // The value is already known, nothing to do.
        } else {
            STORM_LOG_THROW(!normalForm.universalObservationFailureStates.get(initialState), storm::exceptions::NotSupportedException,
                            "Trying to compute undefined conditional probability: the condition has probability 0 under all policies.");
            // The last case for a terminal initial state is that it is already target and the condition is reachable with non-zero probability.
            // In this case, all schedulers induce a conditional probability of 1 (or do not reach the condition, i.e., have undefined value)
            return storm::utility::one<SolutionType>();
        }
    } else {
        // Catch the case where all terminal states have value zero
        if (normalForm.nonZeroTargetStateValues.empty()) {
            return storm::utility::zero<SolutionType>();
        };
    }
    return std::nullopt;  // No trivial case applies, we need to compute the value.
}

}  // namespace internal

template<typename ValueType, typename SolutionType>
std::unique_ptr<CheckResult> computeConditionalProbabilities(Environment const& env, storm::solver::SolveGoal<ValueType, SolutionType>&& goal,
                                                             storm::modelchecker::CheckTask<storm::logic::ConditionalFormula, SolutionType> const& checkTask,
                                                             storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                             storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                             storm::storage::BitVector const& targetStates, storm::storage::BitVector const& conditionStates) {
    // We might require adapting the precision of the solver to counter error propagation (e.g. when computing the normal form).
    auto normalFormConstructionEnv = env;
    auto analysisEnv = env;
    if (env.solver().isForceSoundness()) {
        // We intuitively have to divide the precision into two parts, one for computations when constructing the normal form and one for the actual analysis.
        // As the former is usually less numerically challenging, we use a factor of 1/10 for the normal form construction and 9/10 for the analysis.
        auto const normalFormPrecisionFactor = storm::utility::convertNumber<storm::RationalNumber, std::string>("1/10");
        normalFormConstructionEnv.solver().minMax().setPrecision(env.solver().minMax().getPrecision() * normalFormPrecisionFactor);
        analysisEnv.solver().minMax().setPrecision(env.solver().minMax().getPrecision() *
                                                   (storm::utility::one<storm::RationalNumber>() - normalFormPrecisionFactor));
    }

    // We first translate the problem into a normal form.
    // @see doi.org/10.1007/978-3-642-54862-8_43
    STORM_LOG_THROW(goal.hasRelevantValues(), storm::exceptions::NotSupportedException,
                    "No initial state given. Conditional probabilities can only be computed for models with a single initial state.");
    STORM_LOG_THROW(goal.relevantValues().hasUniqueSetBit(), storm::exceptions::NotSupportedException,
                    "Only one initial state is supported for conditional probabilities");
    STORM_LOG_TRACE("Computing conditional probabilities for a model with " << transitionMatrix.getRowGroupCount() << " states and "
                                                                            << transitionMatrix.getEntryCount() << " transitions.");
    // storm::utility::Stopwatch sw(true);
    auto normalFormData = internal::obtainNormalForm(normalFormConstructionEnv, goal.direction(), transitionMatrix, backwardTransitions, goal.relevantValues(),
                                                     targetStates, conditionStates);
    // sw.stop();
    // STORM_PRINT_AND_LOG("Time for obtaining the normal form: " << sw << ".\n");
    // Then, we solve the induced problem using the selected algorithm
    auto const initialState = *goal.relevantValues().begin();
    ValueType initialStateValue = -storm::utility::one<ValueType>();
    std::unique_ptr<storm::storage::Scheduler<SolutionType>> scheduler = nullptr;
    if (auto trivialValue = internal::handleTrivialCases<ValueType, SolutionType>(initialState, normalFormData); trivialValue.has_value()) {
        initialStateValue = *trivialValue;
        scheduler = std::unique_ptr<storm::storage::Scheduler<SolutionType>>(new storm::storage::Scheduler<SolutionType>(transitionMatrix.getRowGroupCount()));
        STORM_LOG_DEBUG("Initial state has trivial value " << initialStateValue);
    } else {
        STORM_LOG_ASSERT(normalFormData.maybeStates.get(initialState), "Initial state must be a maybe state if it is not a terminal state");
        auto alg = analysisEnv.modelchecker().getConditionalAlgorithmSetting();
        if (alg == ConditionalAlgorithmSetting::Default) {
            alg = ConditionalAlgorithmSetting::Restart;
        }
        STORM_LOG_INFO("Analyzing normal form with " << normalFormData.maybeStates.getNumberOfSetBits() << " maybe states using algorithm '" << alg << ".");
        // sw.restart();
        switch (alg) {
            case ConditionalAlgorithmSetting::Restart: {
                auto result =
                    internal::computeViaRestartMethod(analysisEnv, initialState, goal.direction(), transitionMatrix, backwardTransitions, normalFormData);
                initialStateValue = result.initialStateValue;
                scheduler = std::move(result.scheduler);
                break;
            }
            case ConditionalAlgorithmSetting::Bisection: {
                auto result = internal::computeViaBisection(analysisEnv, internal::BisectionMethodBounds::Simple, initialState, goal, transitionMatrix,
                                                            backwardTransitions, normalFormData);
                initialStateValue = result.initialStateValue;
                scheduler = std::move(result.scheduler);
                break;
            }
            case ConditionalAlgorithmSetting::BisectionAdvanced: {
                auto result = internal::computeViaBisection(analysisEnv, internal::BisectionMethodBounds::Advanced, initialState, goal, transitionMatrix,
                                                            backwardTransitions, normalFormData);
                initialStateValue = result.initialStateValue;
                scheduler = std::move(result.scheduler);
                break;
            }
            case ConditionalAlgorithmSetting::PolicyIteration: {
                initialStateValue = internal::computeViaPolicyIteration(analysisEnv, initialState, goal.direction(), transitionMatrix, normalFormData);
                break;
            }
            default: {
                STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Unknown conditional probability algorithm: " << alg);
            }
        }
        // sw.stop();
        // STORM_PRINT_AND_LOG("Time for analyzing the normal form: " << sw << ".\n");
    }
    std::unique_ptr<CheckResult> result(new ExplicitQuantitativeCheckResult<SolutionType>(initialState, initialStateValue));

    // if produce schedulers was set, we have to construct a scheduler with memory
    if (checkTask.isProduceSchedulersSet() && scheduler) {
        // not sure about this
        storm::utility::graph::computeSchedulerProb1E(normalFormData.targetStates, transitionMatrix, backwardTransitions, normalFormData.targetStates,
                                                      targetStates, *scheduler);
        storm::utility::graph::computeSchedulerProb1E(normalFormData.conditionStates, transitionMatrix, backwardTransitions, normalFormData.conditionStates,
                                                      conditionStates, *scheduler);
        // fill in the scheduler with default choices for states that are missing a choice, these states should be just the ones from which the condition is
        // unreachable this is also used to fill choices for the trivial cases
        for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
            if (!scheduler->isChoiceSelected(state)) {
                // select an arbitrary choice
                scheduler->setChoice(0, state);
            }
        }

        // create scheduler with memory structure
        storm::storage::MemoryStructure::TransitionMatrix memoryTransitions(3, std::vector<boost::optional<storm::storage::BitVector>>(3, boost::none));
        storm::models::sparse::StateLabeling memoryStateLabeling(3);
        memoryStateLabeling.addLabel("init_memory");
        memoryStateLabeling.addLabel("condition_reached");
        memoryStateLabeling.addLabel("target_reached");
        memoryStateLabeling.addLabelToState("init_memory", 0);
        memoryStateLabeling.addLabelToState("condition_reached", 1);
        memoryStateLabeling.addLabelToState("target_reached", 2);

        storm::storage::BitVector allTransitions(transitionMatrix.getEntryCount(), true);
        storm::storage::BitVector conditionExitTransitions(transitionMatrix.getEntryCount(), false);
        storm::storage::BitVector targetExitTransitions(transitionMatrix.getEntryCount(), false);

        for (auto state : conditionStates) {
            for (auto choice : transitionMatrix.getRowGroupIndices(state)) {
                for (auto entryIt = transitionMatrix.getRow(choice).begin(); entryIt < transitionMatrix.getRow(choice).end(); ++entryIt) {
                    conditionExitTransitions.set(entryIt - transitionMatrix.begin(), true);
                }
            }
        }
        for (auto state : targetStates) {
            for (auto choice : transitionMatrix.getRowGroupIndices(state)) {
                for (auto entryIt = transitionMatrix.getRow(choice).begin(); entryIt < transitionMatrix.getRow(choice).end(); ++entryIt) {
                    targetExitTransitions.set(entryIt - transitionMatrix.begin(), true);
                }
            }
        }

        memoryTransitions[0][0] =
            allTransitions & ~conditionExitTransitions & ~targetExitTransitions;  // if neither condition nor target reached, stay in init_memory
        memoryTransitions[0][1] = conditionExitTransitions;
        memoryTransitions[0][2] = targetExitTransitions & ~conditionExitTransitions;
        memoryTransitions[1][1] = allTransitions;  // once condition reached, stay in that memory state
        memoryTransitions[2][2] = allTransitions;  // once target reached, stay in that memory state

        // this assumes there is a single initial state
        auto memoryStructure = storm::storage::MemoryStructure(memoryTransitions, memoryStateLabeling, std::vector<uint64_t>(1, 0), true);

        auto finalScheduler = std::unique_ptr<storm::storage::Scheduler<SolutionType>>(
            new storm::storage::Scheduler<SolutionType>(transitionMatrix.getRowGroupCount(), std::move(memoryStructure)));

        for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
            // set choices for memory 0
            if (conditionStates.get(state)) {
                finalScheduler->setChoice(normalFormData.schedulerChoicesForReachingTargetStates[state], state, 0);
            } else if (targetStates.get(state)) {
                finalScheduler->setChoice(normalFormData.schedulerChoicesForReachingConditionStates[state], state, 0);
            } else {
                finalScheduler->setChoice(scheduler->getChoice(state), state, 0);
            }

            // set choices for memory 1, these are the choices after condition was reached
            finalScheduler->setChoice(normalFormData.schedulerChoicesForReachingTargetStates[state], state, 1);
            // set choices for memory 2, these are the choices after target was reached
            finalScheduler->setChoice(normalFormData.schedulerChoicesForReachingConditionStates[state], state, 2);
        }

        result->asExplicitQuantitativeCheckResult<SolutionType>().setScheduler(std::move(finalScheduler));
    }
    return result;
}

template std::unique_ptr<CheckResult> computeConditionalProbabilities(Environment const& env, storm::solver::SolveGoal<double>&& goal,
                                                                      storm::modelchecker::CheckTask<storm::logic::ConditionalFormula, double> const& checkTask,
                                                                      storm::storage::SparseMatrix<double> const& transitionMatrix,
                                                                      storm::storage::SparseMatrix<double> const& backwardTransitions,
                                                                      storm::storage::BitVector const& targetStates,
                                                                      storm::storage::BitVector const& conditionStates);

template std::unique_ptr<CheckResult> computeConditionalProbabilities(
    Environment const& env, storm::solver::SolveGoal<storm::RationalNumber>&& goal,
    storm::modelchecker::CheckTask<storm::logic::ConditionalFormula, storm::RationalNumber> const& checkTask,
    storm::storage::SparseMatrix<storm::RationalNumber> const& transitionMatrix, storm::storage::SparseMatrix<storm::RationalNumber> const& backwardTransitions,
    storm::storage::BitVector const& targetStates, storm::storage::BitVector const& conditionStates);

}  // namespace storm::modelchecker

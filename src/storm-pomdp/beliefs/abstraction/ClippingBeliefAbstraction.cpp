#include "storm-pomdp/beliefs/abstraction/ClippingBeliefAbstraction.h"

#include "storm-pomdp/beliefs/storage/Belief.h"

#include <optional>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/solver/GlpkLpSolver.h"
#include "storm/solver/LpSolver.h"
#include "storm/storage/expressions/Expression.h"
#include "storm/storage/expressions/ExpressionManager.h"
#include "storm/storage/expressions/Variable.h"
#include "storm/utility/solver.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
ClippingBeliefAbstraction<BeliefType>::ClippingBeliefAbstraction(std::vector<uint64_t>&& observationResolutions)
    : observationResolutions(std::forward<std::vector<uint64_t>>(observationResolutions)) {
    STORM_LOG_ASSERT(std::all_of(observationResolutions.begin(), observationResolutions.end(), [](auto o) { return o > 0; }),
                     "Expected that the resolutions are positive.");
}

template<typename BeliefType>
typename ClippingBeliefAbstraction<BeliefType>::BeliefClipping ClippingBeliefAbstraction<BeliefType>::clipBeliefToGrid(
    const BeliefType& belief, uint64_t resolution, const storm::storage::BitVector& isInfinite) {
    if (!lpSolver) {
        lpSolver = storm::utility::solver::getLpSolver<BeliefValueType>("POMDP LP Solver");
    } else {
        lpSolver->pop();
    }
    lpSolver->push();

    auto const resolutionConverted = storm::utility::convertNumber<BeliefValueType>(resolution);

    std::vector<BeliefValueType> helper(belief.size(), storm::utility::zero<BeliefValueType>());
    helper[0] = resolutionConverted;
    bool done = false;
    // Set-up Variables
    std::vector<storm::expressions::Expression> decisionVariables;
    // Add variable for the clipping value, it is to be minimized
    auto bigDelta = lpSolver->addBoundedContinuousVariable("D", storm::utility::zero<BeliefValueType>(), storm::utility::one<BeliefValueType>(),
                                                           storm::utility::one<BeliefValueType>());
    // State clipping values
    std::vector<storm::expressions::Expression> deltas;
    uint64_t i = 0;
    belief.forEach([this, &isInfinite, &deltas, &i](BeliefStateType const& state, BeliefValueType const& beliefValue) {
        // This is a quite dirty fix to enable GLPK for the TACAS '22 implementation without substantially changing the implementation for Gurobi.
        if (typeid(*lpSolver) == typeid(storm::solver::GlpkLpSolver<BeliefValueType>) && !isInfinite.empty()) {
            if (isInfinite[state]) {
                auto localDelta = lpSolver->addBoundedContinuousVariable("d_" + std::to_string(i), storm::utility::zero<BeliefValueType>(), beliefValue);
                auto deltaExpr = storm::expressions::Expression(localDelta);
                deltas.push_back(deltaExpr);
                lpSolver->addConstraint("state_val_inf_" + std::to_string(i), deltaExpr == lpSolver->getConstant(storm::utility::zero<BeliefValueType>()));
            }
        } else {
            BeliefValueType bound = beliefValue;
            if (!isInfinite.empty()) {
                bound = isInfinite[state] ? storm::utility::zero<BeliefValueType>() : beliefValue;
            }
            auto localDelta = lpSolver->addBoundedContinuousVariable("d_" + std::to_string(i), storm::utility::zero<BeliefValueType>(), bound);
            deltas.push_back(storm::expressions::Expression(localDelta));
        }
        ++i;
    });
    lpSolver->update();
    std::vector<BeliefType> gridCandidates;
    while (!done) {
        BeliefBuilder<BeliefType> candidateBuilder;
        candidateBuilder.setObservation(belief.observation());

        uint64_t j{0};
        uint64_t const jMax = belief.size() - 1;
        belief.forEach([&helper, &j, &resolutionConverted, &candidateBuilder, &jMax](BeliefStateType const& state, BeliefValueType const& beliefValue) {
            if (j < jMax) {
                if (!BeliefNumerics<BeliefValueType>::isZero(helper[j] - helper[j + 1])) {
                    candidateBuilder.addValue(state, (helper[j] - helper[j + 1]) / resolutionConverted);
                }
            } else {
                if (!BeliefNumerics<BeliefValueType>::isZero(helper[jMax])) {
                    candidateBuilder.addValue(state, helper[jMax] / resolutionConverted);
                }
            }
            ++j;
        });
        auto candidate = candidateBuilder.build();
        if (candidate == belief) {
            // TODO Improve handling of successors which are already on the grid
            return BeliefClipping{false, std::move(candidate), storm::utility::zero<BeliefValueType>(), {}, true};
        } else {
            gridCandidates.push_back(candidate);

            // Add variables a_j
            auto decisionVar = lpSolver->addBinaryVariable("a_" + std::to_string(gridCandidates.size() - 1));
            decisionVariables.push_back(storm::expressions::Expression(decisionVar));
            lpSolver->update();

            i = 0;
            belief.forEachCombine(candidate, [&](BeliefStateType const& state, BeliefValueType const& beliefValue, BeliefValueType const& candidateValue) {
                // Add the constraint to describe the transformation between the state values in the beliefs
                // d_i
                storm::expressions::Expression leftSide = deltas[i];
                storm::expressions::Expression targetValue = lpSolver->getConstant(candidateValue);

                // b(s_i) - b_j(s_i) + D * b_j(s_i) - 1 + a_j
                storm::expressions::Expression rightSide =
                    lpSolver->getConstant(beliefValue) - targetValue + storm::expressions::Expression(bigDelta) * targetValue -
                    lpSolver->getConstant(storm::utility::one<BeliefValueType>()) + storm::expressions::Expression(decisionVar);

                // Add left >= right
                lpSolver->addConstraint("state_eq_" + std::to_string(i) + "_" + std::to_string(gridCandidates.size() - 1), leftSide >= rightSide);
                ++i;
                lpSolver->update();
            });
        }
        if (helper.back() == storm::utility::convertNumber<BeliefValueType>(resolution)) {
            // If the last entry of helper is the gridResolution, we have enumerated all necessary distributions
            done = true;
        } else {
            // Update helper by finding the index to increment
            auto helperIt = helper.end() - 1;
            while (*helperIt == *(helperIt - 1)) {
                --helperIt;
            }
            STORM_LOG_ASSERT(helperIt != helper.begin(), "Error in grid clipping - index wrong");
            // Increment the value at the index
            *helperIt += 1;
            // Reset all indices greater than the changed one to 0
            ++helperIt;
            while (helperIt != helper.end()) {
                *helperIt = 0;
                ++helperIt;
            }
        }
    }

    // Only one target belief should be chosen
    lpSolver->addConstraint("choice", storm::expressions::sum(decisionVariables) == lpSolver->getConstant(storm::utility::one<BeliefValueType>()));
    // Link D and d_i
    lpSolver->addConstraint("delta", storm::expressions::Expression(bigDelta) == storm::expressions::sum(deltas));
    // Exclude D = 0 (self-loop)
    lpSolver->addConstraint("not_zero", storm::expressions::Expression(bigDelta) > lpSolver->getConstant(storm::utility::zero<BeliefValueType>()));

    lpSolver->update();

    lpSolver->optimize();
    // Get the optimal belief for clipping
    // Not a belief but has the same type
    BeliefFlatMap<BeliefValueType> deltaValues;
    auto optDelta = storm::utility::zero<BeliefValueType>();
    auto deltaSum = storm::utility::zero<BeliefValueType>();
    if (lpSolver->isOptimal()) {
        uint64_t targetBeliefIndex = std::numeric_limits<uint64_t>::max();
        optDelta = lpSolver->getObjectiveValue();
        for (uint64_t dist = 0; dist < gridCandidates.size(); ++dist) {
            if (lpSolver->getBinaryValue(lpSolver->getManager().getVariable("a_" + std::to_string(dist)))) {
                targetBeliefIndex = dist;
                break;
            }
        }
        STORM_LOG_ASSERT(targetBeliefIndex < gridCandidates.size(), "lp optimal but no belief selected");
        auto targetBelief = gridCandidates.at(targetBeliefIndex);
        i = 0;
        belief.forEachStateInSupport([this, &i, &deltaValues, &deltaSum](BeliefStateType const& state) {
            auto val = lpSolver->getContinuousValue(lpSolver->getManager().getVariable("d_" + std::to_string(i)));
            if (!BeliefNumerics<BeliefValueType>::lessOrEqual(val, storm::utility::zero<BeliefValueType>())) {
                deltaValues.emplace(state, val);
                deltaSum += val;
            }
            ++i;
        });

        if (BeliefNumerics<BeliefValueType>::isZero(optDelta)) {
            // If we get an optimal value of 0, the LP solver considers two beliefs to be equal, possibly due to numerical instability
            // For a sound result, we consider the state to not be clippable
            STORM_LOG_WARN("LP solver returned an optimal value of 0. This should definitely not happen when using a grid");
            STORM_LOG_WARN("Origin" << belief.toString());
            STORM_LOG_WARN("Target [Bel " << targetBelief.toString());
            return BeliefClipping{false, std::move(targetBelief), storm::utility::zero<BeliefValueType>(), {}, false};
        }

        if (optDelta == storm::utility::one<BeliefValueType>()) {
            STORM_LOG_WARN("LP solver returned an optimal value of 1. Sum of state clipping values is " << deltaSum);
            // If we get an optimal value of 1, we cannot clip the belief as by definition this would correspond to a division by 0.
            STORM_LOG_DEBUG("Origin " << belief.toString());
            STORM_LOG_DEBUG("Target " << targetBelief.toString());

            if (deltaSum == storm::utility::one<BeliefValueType>()) {
                return BeliefClipping{false, std::move(targetBelief), storm::utility::zero<BeliefValueType>(), {}, false};
            }
            optDelta = deltaSum;
        }
    }
    return BeliefClipping{lpSolver->isOptimal(), belief, optDelta, deltaValues, false};
}

template class ClippingBeliefAbstraction<Belief<double>>;
template class ClippingBeliefAbstraction<Belief<storm::RationalNumber>>;

}  // namespace storm::pomdp::beliefs
#pragma once
#include <vector>

#include "storm/modelchecker/helper/infinitehorizon/SparseNondeterministicInfiniteHorizonHelper.h"
#include "storm/modelchecker/multiobjective/Objective.h"
#include "storm/modelchecker/multiobjective/pcaa/PcaaWeightVectorChecker.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/solver/LpSolver.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/storage/StronglyConnectedComponentDecomposition.h"
#include "storm/storage/geometry/Halfspace.h"
#include "storm/utility/Extremum.h"
#include "storm/utility/macros.h"
#include "storm/utility/solver.h"
#include "storm/utility/vector.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"

namespace storm {

class Environment;

namespace modelchecker::multiobjective {
template<typename ValueType>
class LongRunAverageSatisfactionHelper {
   public:
    using Point = std::vector<ValueType>;
    using Halfspace = storm::storage::geometry::Halfspace<ValueType>;
    using InfHelperType = storm::modelchecker::helper::SparseNondeterministicInfiniteHorizonHelper<ValueType>;
    LongRunAverageSatisfactionHelper(std::vector<Objective<ValueType>> const& objectives, InfHelperType& helper,
                                     std::vector<std::vector<ValueType>> const& objectiveActionRewards,
                                     std::vector<std::vector<ValueType>> const& objectiveStateRewards,
                                     storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs)
        : helper(helper),
          objectiveActionRewards(objectiveActionRewards),
          objectiveStateRewards(objectiveStateRewards),
          mecs(mecs),
          numInputObjectives(objectives.size()),
          mecAchievableApproximations(mecs.size()) {
        // initialize mec data
        toLocalMecState.assign(helper.getTransitionMatrix().getRowGroupCount(), std::numeric_limits<uint64_t>::max());
        for (auto const& mec : mecs) {
            uint64_t localState = 0;
            for (auto const& [state, _] : mec) {
                toLocalMecState[state] = localState++;
            }
        }

        // initialize data for LRASAT objectives
        minimizingObjectives.resize(numInputObjectives, false);  // will be resized below
        for (uint64_t objIndex = 0; objIndex < objectives.size(); ++objIndex) {
            storm::logic::OperatorFormula const& objFormula = *objectives[objIndex].formula;
            if (objFormula.isProbabilityOperatorFormula() && objFormula.getSubformula().isLongRunAverageRewardFormula()) {
                auto const bound = objFormula.getSubformula().asLongRunAverageRewardFormula().getBound();
                STORM_LOG_THROW(!storm::logic::isStrict(bound.comparisonType), storm::exceptions::NotSupportedException,
                                "Strict bounds are not supported for LRASAT objectives.");
                ValueType threshold = bound.evaluateThresholdAs<ValueType>();
                if (bound.comparisonType == storm::logic::ComparisonType::LessEqual) {
                    minimizingObjectives.set(toInputObjectiveIndexMap.size(), true);  // todo: this currently assumes Pmax=? queries
                    threshold = -threshold;
                }
                satisfactionObjectiveThresholds.push_back(threshold);
                toInputObjectiveIndexMap.push_back(objIndex);
            }
        }
        // initialize data for (expected) LRA objectives
        for (uint64_t objIndex = 0; objIndex < objectives.size(); ++objIndex) {
            storm::logic::OperatorFormula const& objFormula = *objectives[objIndex].formula;
            if (objFormula.isRewardOperatorFormula() && objFormula.getSubformula().isLongRunAverageRewardFormula()) {
                if (storm::solver::minimize(objFormula.getOptimalityType())) {
                    minimizingObjectives.set(toInputObjectiveIndexMap.size(), true);
                }
                toInputObjectiveIndexMap.push_back(objIndex);
            }
        }
        // resize vectors as we now know the number of lra objectives
        uint64_t const numLocalObjectives = toInputObjectiveIndexMap.size();
        minimizingObjectives.resize(numLocalObjectives);
    }

    storm::storage::BitVector getSatisfiableObjectives() const {
        assert(false);
    }

    struct ResultType {
        Point satPoint;
        ValueType optimalWeightedSum;
    };
    ResultType optimize(storm::Environment const& env, std::vector<ValueType> const& inputWeightVector, uint64_t mecIndex, ValueType const& precision) {
        STORM_LOG_ASSERT(inputWeightVector.size() == numInputObjectives, "Weight vector has wrong dimension.");
        STORM_LOG_ASSERT(objectiveActionRewards.size() == numInputObjectives, "Number of action reward vectors does not match number of input objectives.");
        STORM_LOG_ASSERT(objectiveStateRewards.empty() || objectiveStateRewards.size() == numInputObjectives,
                         "Number of state reward vectors does not match number of input objectives.");

        auto& mecApprox = mecAchievableApproximations[mecIndex];
        if (mecApprox.achievablePoints.empty()) {
            initializeMecAchievableApproximation(env, mecIndex);
        }

        std::vector<ValueType> reducedInputWeightVector(toInputObjectiveIndexMap.size());
        storm::utility::vector::selectVectorValues(reducedInputWeightVector, toInputObjectiveIndexMap, inputWeightVector);

        while (true) {
            // Determine the current point and approximation gap
            auto overApprox = getOptimumInPolytope(mecApprox.containingHalfspaces, *mecApprox.maxAbsValue, reducedInputWeightVector);
            auto separator = getSeparatingHalfspace(mecApprox.achievablePoints, overApprox.satPoint);
            //            std::cout << "Current achievable points are:\n";
            //            for (auto const& p : mecApprox.achievablePoints) {
            //                std::cout << "\t" << storm::utility::vector::toString(p) << "\n";
            //            }
            //            std::cout << "Current approximation contains " << mecApprox.containingHalfspaces.size() << " half spaces.\n";
            //            for (auto const& hs : mecApprox.containingHalfspaces) {
            //                std::cout << "\t" << hs.toString() << "\n";
            //            }
            //            std::cout << "Current over-approx point is " << storm::utility::vector::toString(overApprox.satPoint) << " with weighted sum "
            //                      << overApprox.optimalWeightedSum << "\n";
            //            std::cout << "Separating half space is " << separator.toString() << "\n";
            ValueType const gap = separator.distance(overApprox.satPoint);
            STORM_PRINT_AND_LOG("Current gap: " << gap << ", precision: " << precision << "\n");
            if (gap <= precision) {
                // Prepare the result and terminate
                ResultType result;
                result.optimalWeightedSum = overApprox.optimalWeightedSum;
                result.satPoint.assign(numInputObjectives, storm::utility::zero<ValueType>());
                for (uint64_t localObjIndex = 0; localObjIndex < toInputObjectiveIndexMap.size(); ++localObjIndex) {
                    uint64_t const inputObjIndex = toInputObjectiveIndexMap[localObjIndex];
                    if (localObjIndex < satisfactionObjectiveThresholds.size()) {
                        if (overApprox.satPoint[localObjIndex] >= satisfactionObjectiveThresholds[localObjIndex]) {
                            result.satPoint[inputObjIndex] = storm::utility::one<ValueType>();  // todo: -1 for Pmin=? queries ?
                        }
                    } else {
                        result.satPoint[inputObjIndex] = overApprox.satPoint[localObjIndex];  // TODO: calculate with precision
                    }
                }
                //                std::cout << "returning with sat point " << storm::utility::vector::toString(result.satPoint) << "\n";
                return result;
            }
            // refine the approximation
            auto optRes = optimizeExpValues(env, separator.normalVector(), mecIndex);
            for (auto const& val : optRes.satPoint) {
                mecApprox.maxAbsValue &= storm::utility::abs(val);
            }
            mecApprox.achievablePoints.push_back(std::move(optRes.satPoint));
            mecApprox.containingHalfspaces.emplace_back(std::move(separator.normalVector()), std::move(optRes.optimalWeightedSum));
        }
    }

   private:
    ResultType optimizeExpValues(storm::Environment const& env, std::vector<ValueType> const& reducedWeightVector, uint64_t mecIndex) {
        storm::storage::MaximalEndComponent const& mec = mecs[mecIndex];
        RewardGetter rew(*this);
        rew.setWeights(reducedWeightVector);

        // optimize the weighted sum of the objectives
        helper.setOptimizationDirection(storm::solver::OptimizationDirection::Maximize);
        helper.setProduceScheduler(true);
        ResultType result;
        result.optimalWeightedSum = helper.computeLraForComponent(env, rew.weightedStateValueGetter(), rew.weightedActionValueGetter(), mec);

        // std::cout << "Optimal weighted sum in MEC " << mecIndex << " is " << result.optimalWeightedSum << "\n";

        // create the sub-mec that only contains the optimal choices
        storm::storage::MaximalEndComponent subMec;
        {
            storm::storage::SparseMatrixBuilder<ValueType> subMatrixBuilder(mec.size(), mec.size());
            uint64_t subState = 0;
            for (auto const& [state, mecChoices] : mec) {
                uint64_t const choice = helper.getTransitionMatrix().getRowGroupIndices()[state] + helper.getProducedOptimalChoices()[state];
                STORM_LOG_ASSERT(mecChoices.contains(choice), "Optimal choice not part of the MEC.");
                for (auto const& entry : helper.getTransitionMatrix().getRow(choice)) {
                    STORM_LOG_ASSERT(mec.containsState(entry.getColumn()), "MEC transition leads outside of MEC.");
                    subMatrixBuilder.addNextValue(subState, toLocalMecState[entry.getColumn()], entry.getValue());
                }
                ++subState;
            }
            auto subMatrix = subMatrixBuilder.build();
            storm::storage::StronglyConnectedComponentDecompositionOptions options;
            options.areOnlyBottomSccsConsidered = true;
            auto bsccs = storm::storage::StronglyConnectedComponentDecomposition(subMatrix, options);
            STORM_LOG_ASSERT(bsccs.size() == 1, "The sub-MEC induced by the optimal choices has more than one BSCC.");
            storm::storage::StronglyConnectedComponent const& bscc = bsccs[0];
            for (auto const& [state, _] : mec) {
                if (bscc.containsState(toLocalMecState[state])) {
                    uint64_t const choice = helper.getTransitionMatrix().getRowGroupIndices()[state] + helper.getProducedOptimalChoices()[state];
                    subMec.addState(state, {choice});
                }
            }
        }

        // compute the individual objective values in the sub-mec
        helper.setProduceScheduler(false);
        for (uint64_t i = 0; i < reducedWeightVector.size(); ++i) {
            result.satPoint.push_back(helper.computeLraForComponent(env, rew.stateValueGetter(i), rew.actionValueGetter(i), subMec));
        }

        return result;
    }

    void initializeMecAchievableApproximation(storm::Environment const& env, uint64_t const mecIndex) {
        // optimize each objective individually to get an initial set of achievable points
        auto& mecApprox = mecAchievableApproximations[mecIndex];
        for (uint64_t objIndex = 0; objIndex < toInputObjectiveIndexMap.size(); ++objIndex) {
            std::vector<ValueType> weightVector(toInputObjectiveIndexMap.size(), storm::utility::zero<ValueType>());
            weightVector[objIndex] = storm::utility::one<ValueType>();
            auto optRes = optimizeExpValues(env, weightVector, mecIndex);
            for (auto const& val : optRes.satPoint) {
                mecApprox.maxAbsValue &= storm::utility::abs(val);
            }
            mecApprox.achievablePoints.push_back(std::move(optRes.satPoint));
            mecApprox.containingHalfspaces.emplace_back(std::move(weightVector), std::move(optRes.optimalWeightedSum));
        }
    }

    ResultType getOptimumInPolytope(std::vector<Halfspace> const& halfspaces, ValueType const& maxAbsValue, std::vector<ValueType> const& weightVector) const {
        // set up and solve a mixed integer linear program.
        // We use indicator variables (i.e. with domain {0,1}) to encode whether a satisfaction objective is satisfied or not.
        // Setting an indicator variable to 1 yields the weight in the objective but also requires that the sat objective threshold is met.
        auto solver = storm::utility::solver::getLpSolver<ValueType>("");
        solver->setOptimizationDirection(storm::OptimizationDirection::Maximize);
        auto maxValExpr = solver->getConstant(maxAbsValue + maxAbsValue);  // used for indicator constraints encoding
        // add variables and coefficients for optimization function given by the input weight vector
        std::vector<storm::expressions::Variable> valueVariables, indicatorVariables;
        for (uint64_t i = 0; i < weightVector.size(); ++i) {
            if (i < satisfactionObjectiveThresholds.size()) {
                indicatorVariables.push_back(solver->addBinaryVariable("a" + std::to_string(i), weightVector[i]));
                valueVariables.push_back(solver->addUnboundedContinuousVariable("x" + std::to_string(i)));
            } else {
                valueVariables.push_back(solver->addUnboundedContinuousVariable("x" + std::to_string(i), weightVector[i]));
            }
        }
        // add halfspace constraints
        for (auto const& h : halfspaces) {
            solver->addConstraint("", h.toExpression(solver->getManager(), valueVariables));
        }
        // add constraints for satisfaction objectives
        for (uint64_t i = 0; i < satisfactionObjectiveThresholds.size(); ++i) {
            // x_i + max - a_i * max >= threshold_i * a_i
            storm::expressions::Expression lhs = valueVariables[i].getExpression() + maxValExpr - (maxValExpr * indicatorVariables[i].getExpression());
            storm::expressions::Expression rhs = solver->getConstant(satisfactionObjectiveThresholds[i]) * indicatorVariables[i].getExpression();
            solver->addConstraint("", lhs >= rhs);
        }
        solver->update();
        solver->optimize();
        STORM_LOG_THROW(solver->isOptimal(), storm::exceptions::UnexpectedException, "MILP solver did not return an optimal solution.");

        ResultType result;
        for (auto const& var : valueVariables) {
            result.satPoint.push_back(solver->getContinuousValue(var));
        }
        result.optimalWeightedSum = solver->getObjectiveValue();
        return result;
    }

    Halfspace getSeparatingHalfspace(std::vector<Point> const& pointset, Point const& outsidePoint) const {
        // set up a linear program that finds a halfspace that contains all points in pointset but not outsidePoint
        // the distance of outsidePoint to the halfspace is maximized
        // the LP is adapted from 10.1007/978-3-642-33386-6_25
        auto solver = storm::utility::solver::getLpSolver<ValueType>("");
        solver->setOptimizationDirection(storm::OptimizationDirection::Maximize);

        // use distance variable y to encode the largest possible lower bound for  w*outsidePoint - w*p (for all p in the pointset)
        auto distanceVar = solver->addUnboundedContinuousVariable("y", storm::utility::one<ValueType>());
        // add variables for the normal vector w
        std::vector<storm::expressions::Variable> normalVectorVariables;
        std::vector<storm::expressions::Expression> normalVectorVariableExpressions;
        for (uint64_t i = 0; i < outsidePoint.size(); ++i) {
            normalVectorVariables.push_back(
                solver->addBoundedContinuousVariable("w" + std::to_string(i), storm::utility::zero<ValueType>(), storm::utility::one<ValueType>()));
            normalVectorVariableExpressions.push_back(normalVectorVariables.back().getExpression());
        }
        // the normal vector variables shall sum up to 1
        solver->addConstraint("", storm::expressions::sum(normalVectorVariableExpressions) == solver->getConstant(storm::utility::one<ValueType>()));

        // all points must be contained in the halfspace
        for (auto const& p : pointset) {
            std::vector<storm::expressions::Expression> pointConstraint;
            for (uint64_t i = 0; i < outsidePoint.size(); ++i) {
                pointConstraint.push_back(solver->getConstant(outsidePoint[i] - p[i]) * normalVectorVariableExpressions[i]);
            }
            solver->addConstraint("", distanceVar.getExpression() <= storm::expressions::sum(pointConstraint));
        }

        solver->update();
        solver->optimize();
        STORM_LOG_THROW(solver->isOptimal(), storm::exceptions::UnexpectedException, "LP solver did not return an optimal solution.");

        std::vector<ValueType> normalVector;
        for (auto const& var : normalVectorVariables) {
            normalVector.push_back(solver->getContinuousValue(var));
        }
        ValueType const offset = storm::utility::vector::dotProduct(normalVector, outsidePoint) - solver->getContinuousValue(distanceVar);
        return Halfspace(std::move(normalVector), std::move(offset));
    }

    struct RewardGetter {
        using Getter = InfHelperType::ValueGetter;

        RewardGetter(LongRunAverageSatisfactionHelper<ValueType> const& parent)
            : parent(parent), weights(parent.numInputObjectives, storm::utility::zero<ValueType>()) {}

        void setWeights(std::vector<ValueType> const& weightVector) {
            STORM_LOG_ASSERT(weightVector.size() == parent.toInputObjectiveIndexMap.size(), "Weight vector has wrong dimension.");
            for (uint64_t localIndex = 0; localIndex < parent.toInputObjectiveIndexMap.size(); ++localIndex) {
                // in case of minimizing objectives, we actually maximize the negative weighted reward
                ValueType w_i = parent.minimizingObjectives.get(localIndex) ? -weightVector[localIndex] : weightVector[localIndex];
                weights[parent.toInputObjectiveIndexMap[localIndex]] = w_i;
            }
        }

        Getter weightedActionValueGetter() const {
            return [this](uint64_t const& a) {
                auto rew = storm::utility::zero<ValueType>();
                for (auto objIndex : parent.toInputObjectiveIndexMap) {
                    rew += weights[objIndex] * parent.objectiveActionRewards[objIndex][a];
                }
                return rew;
            };
        }

        Getter weightedStateValueGetter() const {
            if (parent.objectiveStateRewards.empty() || std::all_of(parent.toInputObjectiveIndexMap.begin(), parent.toInputObjectiveIndexMap.end(),
                                                                    [this](auto const objIndex) { return parent.objectiveStateRewards[objIndex].empty(); })) {
                return [](uint64_t const&) { return storm::utility::zero<ValueType>(); };
            }
            return [this](uint64_t const& s) {
                auto rew = storm::utility::zero<ValueType>();
                for (auto objIndex : parent.toInputObjectiveIndexMap) {
                    rew += weights[objIndex] * parent.objectiveStateRewards[objIndex][s];
                }
                return rew;
            };
        }

        Getter actionValueGetter(uint64_t const& localObjectiveIndex) {
            uint64_t const liftedObjIndex = parent.toInputObjectiveIndexMap[localObjectiveIndex];
            auto const& rewVec = parent.objectiveActionRewards[liftedObjIndex];
            if (parent.minimizingObjectives.get(localObjectiveIndex)) {
                return [&rewVec](uint64_t const& a) { return -rewVec[a]; };

            } else {
                return [&rewVec](uint64_t const& a) { return rewVec[a]; };
            }
        }

        Getter stateValueGetter(uint64_t const& localObjectiveIndex) {
            uint64_t const liftedObjIndex = parent.toInputObjectiveIndexMap[localObjectiveIndex];
            if (parent.objectiveStateRewards.empty() || parent.objectiveStateRewards[liftedObjIndex].empty()) {
                return [](uint64_t const&) { return storm::utility::zero<ValueType>(); };
            }
            auto const& rewVec = parent.objectiveStateRewards[liftedObjIndex];
            if (parent.minimizingObjectives.get(localObjectiveIndex)) {
                return [&rewVec](uint64_t const& s) { return -rewVec[s]; };
            } else {
                return [&rewVec](uint64_t const& s) { return rewVec[s]; };
            }
        }

        LongRunAverageSatisfactionHelper<ValueType> const& parent;
        std::vector<ValueType> weights;
    };

    struct AchApproximation {
        std::vector<Point> achievablePoints;             // all points are achievable
        std::vector<Halfspace> containingHalfspaces;     // all halfspaces contain the achievable points
        storm::utility::Maximum<ValueType> maxAbsValue;  // upper bound on the maximum absolute value of any objective in this MEC (used for MILP encoding)
    };

    InfHelperType& helper;
    std::vector<std::vector<ValueType>> const& objectiveActionRewards;
    std::vector<std::vector<ValueType>> const& objectiveStateRewards;
    storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs;
    std::vector<uint64_t> toLocalMecState;  // maps global state indices to local MEC state indices

    uint64_t const numInputObjectives;
    std::vector<uint64_t> toInputObjectiveIndexMap;  // maps local objective indices to input objective indices

    storm::storage::BitVector minimizingObjectives;
    std::vector<ValueType> satisfactionObjectiveThresholds;

    std::vector<AchApproximation> mecAchievableApproximations;  // indexed by MEC index
};

}  // namespace modelchecker::multiobjective
}  // namespace storm

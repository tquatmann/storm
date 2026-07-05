#pragma once

#include <random>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/exceptions/InvalidOperationException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/transformer/ContinuousToDiscreteTimeModelTransformer.h"
#include "storm/transformer/NonMarkovianChainTransformer.h"
#include "storm/transformer/StatePermuter.h"
#include "storm/transformer/SymbolicToSparseTransformer.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/builder.h"
#include "storm/utility/macros.h"
#include "storm/utility/permutation.h"

namespace storm {
namespace api {

/*!
 * Eliminates chains of non-Markovian states from a given Markov Automaton
 */
template<typename ValueType>
std::pair<std::shared_ptr<storm::models::sparse::Model<ValueType>>, std::vector<std::shared_ptr<storm::logic::Formula const>>> eliminateNonMarkovianChains(
    std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> const& ma, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
    storm::transformer::EliminationLabelBehavior labelBehavior) {
    auto newFormulas = storm::transformer::NonMarkovianChainTransformer<ValueType>::checkAndTransformFormulas(formulas);
    STORM_LOG_WARN_COND(newFormulas.size() == formulas.size(), "The state elimination does not preserve all properties.");
    STORM_LOG_WARN_COND(!(labelBehavior == storm::transformer::EliminationLabelBehavior::KeepLabels),
                        "Labels are not preserved by the state elimination. This may cause incorrect results.");
    return std::make_pair(storm::transformer::NonMarkovianChainTransformer<ValueType>::eliminateNonmarkovianStates(ma, labelBehavior), newFormulas);
}

/*!
 * Transforms the given continuous model to a discrete time model.
 * If such a transformation does not preserve one of the given formulas, a warning is issued.
 */
template<typename ValueType>
std::pair<std::shared_ptr<storm::models::sparse::Model<ValueType>>, std::vector<std::shared_ptr<storm::logic::Formula const>>>
transformContinuousToDiscreteTimeSparseModel(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                             std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) {
    storm::transformer::ContinuousToDiscreteTimeModelTransformer<ValueType> transformer;

    std::string timeRewardName = "_time";
    while (model->hasRewardModel(timeRewardName)) {
        timeRewardName += "_";
    }
    auto newFormulas = transformer.checkAndTransformFormulas(formulas, timeRewardName);
    STORM_LOG_WARN_COND(newFormulas.size() == formulas.size(),
                        "Transformation of a " << model->getType() << " to a discrete time model does not preserve all properties.");

    if (model->isOfType(storm::models::ModelType::Ctmc)) {
        return std::make_pair(transformer.transform(*model->template as<storm::models::sparse::Ctmc<ValueType>>(), timeRewardName), newFormulas);
    } else if (model->isOfType(storm::models::ModelType::MarkovAutomaton)) {
        return std::make_pair(transformer.transform(*model->template as<storm::models::sparse::MarkovAutomaton<ValueType>>(), timeRewardName), newFormulas);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                        "Transformation of a " << model->getType() << " to a discrete time model is not supported");
    }
    return std::make_pair(nullptr, newFormulas);
}

/*!
 * Transforms the given continuous model to a discrete time model IN PLACE.
 * This means that the input continuous time model is replaced by the new discrete time model.
 * If such a transformation does not preserve one of the given formulas, an error is issued.
 */
template<typename ValueType>
std::pair<std::shared_ptr<storm::models::sparse::Model<ValueType>>, std::vector<std::shared_ptr<storm::logic::Formula const>>>
transformContinuousToDiscreteTimeSparseModel(storm::models::sparse::Model<ValueType>&& model,
                                             std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) {
    storm::transformer::ContinuousToDiscreteTimeModelTransformer<ValueType> transformer;

    std::string timeRewardName = "_time";
    while (model.hasRewardModel(timeRewardName)) {
        timeRewardName += "_";
    }
    auto newFormulas = transformer.checkAndTransformFormulas(formulas, timeRewardName);
    STORM_LOG_WARN_COND(newFormulas.size() == formulas.size(),
                        "Transformation of a " << model.getType() << " to a discrete time model does not preserve all properties.");

    if (model.isOfType(storm::models::ModelType::Ctmc)) {
        return std::make_pair(transformer.transform(std::move(*model.template as<storm::models::sparse::Ctmc<ValueType>>()), timeRewardName), newFormulas);
    } else if (model.isOfType(storm::models::ModelType::MarkovAutomaton)) {
        return std::make_pair(transformer.transform(std::move(*model.template as<storm::models::sparse::MarkovAutomaton<ValueType>>()), timeRewardName),
                              newFormulas);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                        "Transformation of a " << model.getType() << " to a discrete time model is not supported.");
    }
    return std::make_pair(nullptr, newFormulas);
}

template<typename ValueType>
std::shared_ptr<storm::logic::Formula const> checkAndTransformContinuousToDiscreteTimeFormula(storm::logic::Formula const& formula,
                                                                                              std::string const& timeRewardName = "_time") {
    storm::transformer::ContinuousToDiscreteTimeModelTransformer<ValueType> transformer;
    if (transformer.preservesFormula(formula)) {
        return transformer.checkAndTransformFormulas({formula.asSharedPointer()}, timeRewardName).front();
    } else {
        STORM_LOG_ERROR("Unable to transform formula " << formula << " to discrete time.");
    }
    return nullptr;
}

/*!
 * Transforms the given symbolic model to a sparse model.
 */
template<storm::dd::DdType Type, typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> transformSymbolicToSparseModel(
    std::shared_ptr<storm::models::symbolic::Model<Type, ValueType>> const& symbolicModel,
    std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas = std::vector<std::shared_ptr<storm::logic::Formula const>>()) {
    switch (symbolicModel->getType()) {
        case storm::models::ModelType::Dtmc:
            return storm::transformer::SymbolicDtmcToSparseDtmcTransformer<Type, ValueType>().translate(
                *symbolicModel->template as<storm::models::symbolic::Dtmc<Type, ValueType>>(), formulas);
        case storm::models::ModelType::Mdp:
            return storm::transformer::SymbolicMdpToSparseMdpTransformer<Type, ValueType>::translate(
                *symbolicModel->template as<storm::models::symbolic::Mdp<Type, ValueType>>(), formulas);
        case storm::models::ModelType::Ctmc:
            return storm::transformer::SymbolicCtmcToSparseCtmcTransformer<Type, ValueType>::translate(
                *symbolicModel->template as<storm::models::symbolic::Ctmc<Type, ValueType>>(), formulas);
        case storm::models::ModelType::MarkovAutomaton:
            return storm::transformer::SymbolicMaToSparseMaTransformer<Type, ValueType>::translate(
                *symbolicModel->template as<storm::models::symbolic::MarkovAutomaton<Type, ValueType>>(), formulas);
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                            "Transformation of symbolic " << symbolicModel->getType() << " to sparse model is not supported.");
    }
    return nullptr;
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> transformToNondeterministicModel(storm::models::sparse::Model<ValueType>&& model) {
    storm::storage::sparse::ModelComponents<ValueType> components(std::move(model.getTransitionMatrix()), std::move(model.getStateLabeling()),
                                                                  std::move(model.getRewardModels()));
    components.choiceLabeling = std::move(model.getOptionalChoiceLabeling());
    components.stateValuations = std::move(model.getOptionalStateValuations());
    components.choiceOrigins = std::move(model.getOptionalChoiceOrigins());
    if (model.isOfType(storm::models::ModelType::Dtmc)) {
        return storm::utility::builder::buildModelFromComponents(storm::models::ModelType::Mdp, std::move(components));
    } else if (model.isOfType(storm::models::ModelType::Ctmc)) {
        components.rateTransitions = true;
        components.markovianStates = storm::storage::BitVector(components.transitionMatrix.getRowGroupCount(), true);
        return storm::utility::builder::buildModelFromComponents(storm::models::ModelType::MarkovAutomaton, std::move(components));
    } else {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException,
                        "Cannot transform model of type " << model.getType() << " to a nondeterministic model.");
    }
}

/*!
 * Permutes the order of the states of the model according to the given order.
 * The order of the available choices at a state (of a nondeterministic model) is not changed.
 * A seed can be given which will be respected if a random permutation is requested.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> permuteModelStates(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                                                            storm::utility::permutation::OrderKind order,
                                                                            std::optional<uint64_t> seed = std::nullopt) {
    std::vector<storm::utility::permutation::index_type> permutation;
    if (order == storm::utility::permutation::OrderKind::Random && seed.has_value()) {
        permutation = storm::utility::permutation::createRandomPermutation(model->getNumberOfStates(), seed.value());
    } else {
        permutation = storm::utility::permutation::createPermutation(order, model->getTransitionMatrix(), model->getInitialStates());
    }
    return storm::transformer::permuteStates(*model, permutation);
}

/*!
 * Perturbs the outgoing transitions of the given model with '1 - gamma' probability up to 'delta' in L_1 distance, and with probability 'gamma' up to '2 *
 * delta', by shifting the intervals. This is inspired by the perturbation that Kiefer & Tang apply in their paper 'Approximate Bisimulation Minimization':
 * https://doi.org/10.48550/arXiv.2110.00326
 * Note that this is an in-place operation, i.e., the perturbation gets applied directly on the input model.
 * @param model input model, must contain interval uncertainty.
 * @param delta amount of perturbation per (choice/action of) state given as L_1 distance over the outgoing distribution.
 * @param gamma probability to perturb the outgoing transitions of a state up to '2 * delta'.
 * @return perturbed input model.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> perturbModelTransitions(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                                                                 double delta, double gamma) {
    using SolutionType = storm::IntervalBaseType<ValueType>;
    auto const zero = storm::utility::zero<SolutionType>();
    auto const one = storm::utility::one<SolutionType>();

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> pickUniformBetween01(0.0, 1.0);
    std::bernoulli_distribution pickBernoulliForGamma(gamma);

    auto numberOfPerturbedRows = 0;
    auto stateChoiceIndices = model->getTransitionMatrix().getRowGroupIndices();

    for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
        // Perturb every possible choice/action of a state.
        auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
        for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
            auto currentChoice = stateChoiceIndices[s] + choiceOffset;
            auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);
            const auto outDegreeOfCurrentChoice = static_cast<std::size_t>(end - begin);

            // Cannot perturb transitions if there are less than two.
            if (outDegreeOfCurrentChoice < 2) {
                // Nothing to do here.
                continue;
            }

            // Choose the amount of perturbation with respect to the error rate.
            SolutionType remainingL1 = pickBernoulliForGamma(rng) ? 2.0 * delta : delta;

            // Choose how often to try before admitting row is too tight.
            int retries = 2000;
            std::uniform_int_distribution<int> pickIndex(0, static_cast<int>(outDegreeOfCurrentChoice) - 1);

            if constexpr (!storm::IsIntervalType<ValueType>) {
                // Restructure transitions in contagious data structures.
                std::vector<SolutionType> probabilities;
                probabilities.reserve(outDegreeOfCurrentChoice);

                for (auto it = begin; it != end; ++it) {
                    probabilities.push_back(it->getValue());
                }

                while (remainingL1 > 1e-15 && retries-- > 0) {
                    // Pick two different target states.
                    int i = pickIndex(rng), j = pickIndex(rng);
                    if (i == j)
                        continue;

                    // Move mass theta from transition to target state 'j' to transition of target state 'i'.
                    // Make sure that we do not shift more mass than P('s', 'currentChoice', 'j') provides and do not exceed the headroom that P('s',
                    // 'currentChoice', 'i') offers.
                    SolutionType headroomAtI = 1 - probabilities[i];
                    SolutionType cap = std::min(probabilities[j], headroomAtI);
                    if (cap <= zero)
                        continue;

                    // Spend at most remaining L_1/2 (because this step costs 2 * theta).
                    SolutionType halfRemainingL1 = remainingL1 / SolutionType(2);
                    SolutionType theta = std::min(cap, halfRemainingL1);

                    // Randomize within [0, theta].
                    theta *= pickUniformBetween01(rng);

                    if (theta > zero) {
                        probabilities[i] += theta;
                        probabilities[j] -= theta;
                        remainingL1 = remainingL1 - (2.0 * theta);
                    }
                }

                auto it = begin;
                for (std::size_t k = 0; k < outDegreeOfCurrentChoice; ++k, ++it) {
                    auto p = std::max(zero, std::min(one, probabilities[k]));
                    it->setValue(p);
                }

                numberOfPerturbedRows++;
            } else {
                // Restructure lower and upper bounds of transitions in contagious data structures.
                std::vector<SolutionType> lowerBounds;
                std::vector<SolutionType> upperBounds;
                lowerBounds.reserve(outDegreeOfCurrentChoice);
                upperBounds.reserve(outDegreeOfCurrentChoice);

                for (auto it = begin; it != end; ++it) {
                    lowerBounds.push_back(it->getValue().lower());
                    upperBounds.push_back(it->getValue().upper());
                }

                while (remainingL1 > 1e-15 && retries-- > 0) {
                    // Pick two different states.
                    int i = pickIndex(rng), j = pickIndex(rng);
                    if (i == j)
                        continue;

                    // Move mass theta from state (choice) j -> i on both lowerBounds and upperBounds.
                    // Add at i (both <= 1), subtract at j (both >= 0).
                    SolutionType capPlus = std::min(one - upperBounds[i], one - lowerBounds[i]);
                    SolutionType capMinus = std::min(upperBounds[j], lowerBounds[j]);
                    SolutionType thetaMax = std::min(capPlus, capMinus);
                    if (thetaMax <= zero)
                        continue;

                    // Spend at most remaining L_1/2 (because this step costs 2 * theta).
                    SolutionType halfRemainingL1 = remainingL1 / SolutionType(2);
                    SolutionType theta = std::min(thetaMax, halfRemainingL1);
                    // Randomize within [0, theta].
                    theta *= pickUniformBetween01(rng);

                    if (theta > zero) {
                        lowerBounds[i] += theta;
                        upperBounds[i] += theta;
                        lowerBounds[j] -= theta;
                        upperBounds[j] -= theta;
                        remainingL1 = remainingL1 - (2.0 * theta);
                    }
                }

                // Clamp tiny numeric drifts and write back to model.
                auto it = begin;
                for (std::size_t k = 0; k < outDegreeOfCurrentChoice; ++k, ++it) {
                    auto l = std::max(zero, std::min(one, lowerBounds[k]));
                    auto u = std::max(l, std::min(one, upperBounds[k]));

                    STORM_LOG_THROW(l != zero, storm::exceptions::InvalidArgumentException,
                                    "Perturbation of interval " << it->getValue() << " yielded lower bound of zero and thus violates graph-preservation.");

                    it->setValue(ValueType(l, u));
                }

                numberOfPerturbedRows++;
            }
        }
    }

    STORM_PRINT_AND_LOG("Perturbed " << numberOfPerturbedRows << " rows in transition matrix.\n");

    return model;
}

/*!
 * Computes for every transition interval in the model its interval of feasible values, i.e., the values that can actually be instantiated w.r.t. other
 * intervals under the simplex constraint.
 * @param model input model, must be an interval model.
 * @return input model with only feasible transition intervals.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> makeTransitionIntervalsFeasible(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                                                                         double precision) {
    if constexpr (storm::IsIntervalType<ValueType>) {
        using SolutionType = storm::IntervalBaseType<ValueType>;
        auto const zero = storm::utility::zero<SolutionType>();
        auto const one = storm::utility::one<SolutionType>();

        storm::utility::ConstantsComparator<SolutionType> comparator(precision);

        auto numberOfAffectedRows = 0;
        auto stateChoiceIndices = model->getTransitionMatrix().getRowGroupIndices();

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            // Make transition intervals feasible for every possible choice/action of a state.
            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);

                // First, we sum up all the transition intervals in \bar{I}.
                SolutionType lowerSum = zero;
                SolutionType upperSum = zero;

                for (auto it = begin; it != end; ++it) {
                    auto const& currentInterval = it->getValue();
                    lowerSum += currentInterval.lower();
                    upperSum += currentInterval.upper();
                }

                // Feasibility check for the row.
                STORM_LOG_THROW(!(comparator.isLess(one, lowerSum) || comparator.isLess(upperSum, one)), storm::exceptions::InvalidArgumentException,
                                "Choice row " << currentChoice << " has no feasible distribution.");

                // Next, we compute the actual feasible intervals.
                for (auto it = begin; it != end; ++it) {
                    ValueType intervalToCurrentSuccessor = it->getValue();  // Interval I
                    SolutionType lowerOtherStates = lowerSum - intervalToCurrentSuccessor.lower();
                    SolutionType upperOtherStates = upperSum - intervalToCurrentSuccessor.upper();

                    // Represents the possible interval when considering the lower and upper bounds of the other states in the (choice) row:
                    // [1 - \sup (\bar{I} - I), 1 - \inf (\bar{I} - I)]
                    SolutionType lowerImpliedByOtherStates = one - upperOtherStates;
                    SolutionType upperImpliedByOtherStates = one - lowerOtherStates;
                    SolutionType feasibleLower = std::max(intervalToCurrentSuccessor.lower(), lowerImpliedByOtherStates);
                    SolutionType feasibleUpper = std::min(intervalToCurrentSuccessor.upper(), upperImpliedByOtherStates);

                    // [1 - \sup (\bar{I} - I), 1 - \inf (\bar{I} - I)] \cap [0, 1]
                    feasibleLower = std::max(feasibleLower, zero);
                    feasibleUpper = std::min(feasibleUpper, one);

                    // In case the feasible interval should be a point-interval we have to ensure that it is non-empty in the non-exact case.
                    if (comparator.isEqual(feasibleLower, feasibleUpper)) {
                        feasibleUpper = feasibleLower;
                    }

                    ValueType feasibleIntervalToCurrentSuccessor(feasibleLower, carl::BoundType::WEAK, feasibleUpper, carl::BoundType::WEAK);

                    STORM_LOG_THROW(!feasibleIntervalToCurrentSuccessor.isEmpty(), storm::exceptions::InvalidArgumentException,
                                    "Interval " << intervalToCurrentSuccessor << " of choice/state " << currentChoice << " to state " << it->getColumn()
                                                << " yields empty feasible interval.");

                    STORM_LOG_THROW(!comparator.isEqual(feasibleIntervalToCurrentSuccessor.lower(), zero), storm::exceptions::InvalidArgumentException,
                                    "Interval " << intervalToCurrentSuccessor << " of choice/state " << currentChoice << " to state " << it->getColumn()
                                                << " yields non-graph-preserving feasible interval.");

                    // Only replace intervals if the values really changed and are not only due to imprecision.
                    if (!comparator.isEqual(intervalToCurrentSuccessor.lower(), feasibleIntervalToCurrentSuccessor.lower()) ||
                        !comparator.isEqual(intervalToCurrentSuccessor.upper(), feasibleIntervalToCurrentSuccessor.upper())) {
                        it->setValue(feasibleIntervalToCurrentSuccessor);
                        numberOfAffectedRows++;
                    }
                }
            }
        }

        STORM_PRINT_AND_LOG("Computing feasible transition intervals affected " << numberOfAffectedRows << " rows in transition matrix.\n");
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Cannot compute feasible intervals on non-interval model.");
    }

    return model;
}

/*!
 * Creates a point-interval model variant of the input model, i.e., turns each transition probability p into an interval [p, p].
 * TODO: Implement copying of RewardModel.
 * @param model input model, must be a non-interval model.
 * @return input model with transition point-intervals.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<storm::Interval>> transformToPointIntervalModel(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model) {
    if constexpr (storm::IsIntervalType<ValueType>) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Cannot convert interval model to point-interval model.");
    } else {
        if (!model->isOfType(storm::models::ModelType::Dtmc) && !model->isOfType(storm::models::ModelType::Mdp)) {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "We only support converting to point-interval for DTMCs and MDPs.");
        }

        using IntervalType = storm::Interval;
        using RewardModelType = storm::models::sparse::StandardRewardModel<IntervalType>;

        auto const& oldMatrix = model->getTransitionMatrix();
        auto const& rowGroupIndices = oldMatrix.getRowGroupIndices();

        bool const isDeterministicModel = model->isOfType(storm::models::ModelType::Dtmc);
        storm::storage::SparseMatrixBuilder<IntervalType> builder(oldMatrix.getRowCount(), oldMatrix.getColumnCount(), oldMatrix.getNonzeroEntryCount(), true,
                                                                  !isDeterministicModel, isDeterministicModel ? 0 : oldMatrix.getRowGroupCount());

        for (std::size_t s = 0; s < model->getNumberOfStates(); s++) {
            if (!isDeterministicModel) {
                builder.newRowGroup(rowGroupIndices[s]);
            }

            for (auto row = rowGroupIndices[s]; row < rowGroupIndices[s + 1]; row++) {
                for (auto const& entry : oldMatrix.getRow(row)) {
                    ValueType p = storm::utility::convertNumber<ValueType>(entry.getValue());
                    builder.addNextValue(row, entry.getColumn(), IntervalType(p, p));
                }
            }
        }

        storm::storage::sparse::ModelComponents<IntervalType, RewardModelType> components(builder.build(),
                                                                                          storm::models::sparse::StateLabeling(model->getStateLabeling()));

        if (model->hasChoiceLabeling()) {
            components.choiceLabeling = storm::models::sparse::ChoiceLabeling(model->getChoiceLabeling());
        }

        if (model->hasStateValuations()) {
            components.stateValuations = model->getStateValuations();
        }

        if (model->hasChoiceOrigins()) {
            components.choiceOrigins = model->getChoiceOrigins();
        }

        // TODO: Implement copying of RewardModel.

        if (isDeterministicModel) {
            return std::make_shared<storm::models::sparse::Dtmc<IntervalType>>(std::move(components));
        }

        return std::make_shared<storm::models::sparse::Mdp<IntervalType>>(std::move(components));
    }
}

}  // namespace api
}  // namespace storm

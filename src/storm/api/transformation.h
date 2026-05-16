#pragma once

#include <random>

#include "storm/exceptions/InvalidOperationException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/transformer/ContinuousToDiscreteTimeModelTransformer.h"
#include "storm/transformer/NonMarkovianChainTransformer.h"
#include "storm/transformer/StatePermuter.h"
#include "storm/transformer/SymbolicToSparseTransformer.h"
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
 * delta'. Inspired by the perturbation that Kiefer & Tang apply in their paper 'Approximate Bisimulation Minimization':
 * https://doi.org/10.48550/arXiv.2110.00326
 * @param model input model, must contain interval uncertainty.
 * @param delta amount of perturbation per state given as L_1 distance over the outgoing distribution.
 * @param gamma probability to perturb the outgoing transitions of a state up to '2 * delta'.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> perturbModelTransitions(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                                                                 double delta, double gamma) {
    if constexpr (!storm::IsIntervalType<ValueType>) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "No support for perturbation of model for non-interval models yet.");
    }

    using SolutionType = storm::IntervalBaseType<ValueType>;
    auto& transMatrix = const_cast<storm::storage::SparseMatrix<ValueType>&>(model->getTransitionMatrix());

    auto rowCount = transMatrix.getRowCount();

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> pickBetween01(0.0, 1.0);
    std::bernoulli_distribution pickErrorRate(gamma);

    for (std::size_t s = 0; s < rowCount; ++s) {
        auto begin = transMatrix.begin(s), end = transMatrix.end(s);
        const auto outdeg = static_cast<std::size_t>(end - begin);
        if (outdeg == 0)
            continue;

        // Restructure lower and upper bounds of transitions in contagious data structures.
        std::vector<SolutionType> lowerBounds;
        std::vector<SolutionType> upperBounds;
        lowerBounds.reserve(outdeg);
        upperBounds.reserve(outdeg);

        for (auto it = begin; it != end; ++it) {
            lowerBounds.push_back(it->getValue().lower());
            upperBounds.push_back(it->getValue().upper());
        }

        // Choose the amount of perturbation with respect to the error rate.
        double remainingL1 = pickErrorRate(rng) ? 2.0 * delta : delta;

        // Cannot perturb transitions if there are less than two.
        if (outdeg < 2) {
            continue;
        }

        // Choose how often to try before admitting row is too tight.
        int retries = 2000;
        std::uniform_int_distribution<int> pickIndex(0, static_cast<int>(outdeg) - 1);

        while (remainingL1 > 1e-15 && retries-- > 0) {
            // Pick two different states.
            int i = pickIndex(rng), j = pickIndex(rng);
            if (i == j)
                continue;

            // Move mass theta from state j -> i on both lowerBounds and upperBounds.
            // Add at i (both <= 1), subtract at j (both >= 0).
            SolutionType capPlus = std::min(1.0 - upperBounds[i], 1.0 - lowerBounds[i]);
            SolutionType capMinus = std::min(upperBounds[j], lowerBounds[j]);
            SolutionType thetaMax = std::min(capPlus, capMinus);
            if (thetaMax <= 0.0)
                continue;

            // Spend at most remaining L_1/2 (because this step costs 2 * theta).
            double theta = std::min(thetaMax, remainingL1 / 2.0);
            // Randomize within [0, theta].
            theta *= pickBetween01(rng);

            if (theta > 0.0) {
                lowerBounds[i] += theta;
                upperBounds[i] += theta;
                lowerBounds[j] -= theta;
                upperBounds[j] -= theta;
                remainingL1 -= 2.0 * theta;
            }
        }

        // Clamp tiny numeric drifts and write back to model.
        auto it = begin;
        for (std::size_t k = 0; k < outdeg; ++k, ++it) {
            auto l = std::max(0.0, std::min(1.0, lowerBounds[k]));
            auto u = std::max(l, std::min(1.0, upperBounds[k]));
            it->setValue(ValueType(l, u));
        }
    }

    return model;
}

}  // namespace api
}  // namespace storm

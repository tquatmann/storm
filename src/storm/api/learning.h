#pragma once

#include <random>

#include <boost/math/distributions/beta.hpp>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/Distribution.h"

namespace storm {
namespace api {

template<typename ValueType>
storm::Interval getClopperPearsonInterval(std::uint_fast64_t k, std::uint_fast64_t n, double delta) {
    ValueType lowerBound;
    ValueType upperBound;

    if (n < 1) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Tried to obtain clopper-pearson interval for state-action pair with no samples.");
    }

    if (k > n) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException,
                        "Tried to obtain clopper-pearsin interval for state-action pair with more successes than samples. How?");
    }

    if (k > 0) {
        lowerBound = boost::math::ibeta_inv(k, n - k + 1, delta / 2.0);
    } else {
        lowerBound = storm::utility::zero<ValueType>();
    }

    if (k < n) {
        upperBound = boost::math::ibetac_inv(k + 1, n - k, delta / 2.0);
    } else {
        upperBound = storm::utility::one<ValueType>();
    }

    return storm::Interval(lowerBound, upperBound);
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<storm::Interval>> learnIMDPFromMDPByClopperPearson(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, double delta, std::size_t numberOfSamples) {
    if constexpr (!storm::IsIntervalType<ValueType>) {
        using IntervalType = storm::Interval;
        using RewardModelType = storm::models::sparse::StandardRewardModel<IntervalType>;

        std::unordered_map<std::size_t, std::uint_fast64_t> choiceSampleCount;
        std::unordered_map<std::size_t, std::unordered_map<storm::storage::sparse::state_type, std::uint_fast64_t>> choiceSuccessorSampleCount;

        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

        auto stateChoiceIndices = model->getTransitionMatrix().getRowGroupIndices();

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            // Make sure to sample each state-action pair up to 'numberOfSamples'-times.
            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);

                // Construct distribution over successor state of current state-action pair.
                storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> successorDistribution;
                auto successorState = begin;
                while (successorState != end) {
                    successorDistribution.addProbability(successorState->getColumn(), successorState->getValue());
                    successorState++;
                }

                // Draw samples from successor distribution.
                while (choiceSampleCount[currentChoice] < numberOfSamples) {
                    auto sampledSuccessorState = successorDistribution.sampleFromDistribution(uniformDistribution(rng));

                    // Update sample counts.
                    choiceSuccessorSampleCount[currentChoice][sampledSuccessorState] = choiceSuccessorSampleCount[currentChoice][sampledSuccessorState] + 1;
                    choiceSampleCount[currentChoice] = choiceSampleCount[currentChoice] + 1;
                }
            }
        }

        if (!model->isOfType(storm::models::ModelType::Dtmc) && !model->isOfType(storm::models::ModelType::Mdp)) {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "We only support converting to point-interval for DTMCs and MDPs.");
        }

        auto const& oldMatrix = model->getTransitionMatrix();

        storm::storage::SparseMatrixBuilder<IntervalType> builder(oldMatrix.getRowCount(), oldMatrix.getColumnCount(), oldMatrix.getNonzeroEntryCount(), true,
                                                                  true, oldMatrix.getRowGroupCount());

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            builder.newRowGroup(stateChoiceIndices[s]);

            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);

                auto successorState = begin;
                while (successorState != end) {
                    // Compute Clopper-Pearson intervals.
                    IntervalType estimatedInterval = storm::api::getClopperPearsonInterval<ValueType>(
                        choiceSuccessorSampleCount[currentChoice][successorState->getColumn()], choiceSampleCount[currentChoice], delta);

                    builder.addNextValue(currentChoice, successorState->getColumn(), estimatedInterval);

                    successorState++;
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

        if (model->isOfType(storm::models::ModelType::Dtmc)) {
            return std::make_shared<storm::models::sparse::Dtmc<IntervalType>>(std::move(components));
        }

        return std::make_shared<storm::models::sparse::Mdp<IntervalType>>(std::move(components));
    }
}

}  // namespace api
}  // namespace storm
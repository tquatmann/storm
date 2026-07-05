#pragma once

#include <random>

#include <boost/math/distributions/beta.hpp>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/Distribution.h"

namespace storm {
namespace api {

template<typename Iterator>
std::size_t countSuccessors(Iterator begin, Iterator end) {
    std::size_t result = 0;
    while (begin != end) {
        ++result;
        ++begin;
    }
    return result;
}

inline double getComponentFailureProbability(double distributionFailureProbability, std::size_t numberOfSuccessors) {
    STORM_LOG_THROW(numberOfSuccessors > 0, storm::exceptions::InvalidOperationException, "Tried to learn a choice without successors.");
    return distributionFailureProbability / static_cast<double>(numberOfSuccessors);
}

template<typename ValueType>
storm::Interval getClopperPearsonInterval(std::uint_fast64_t k, std::uint_fast64_t n, double componentFailureProbability) {
    ValueType lowerBound;
    ValueType upperBound;

    if (n < 1) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Tried to obtain Clopper-Pearson interval for state-action pair with no samples.");
    }

    if (k > n) {
        STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException,
                        "Tried to obtain Clopper-Pearson interval for state-action pair with more successes than samples. How?");
    }

    if (k > 0) {
        lowerBound = boost::math::ibeta_inv(k, n - k + 1, componentFailureProbability / 2.0);
    } else {
        lowerBound = storm::utility::zero<ValueType>();
    }

    if (k < n) {
        upperBound = boost::math::ibetac_inv(k + 1, n - k, componentFailureProbability / 2.0);
    } else {
        upperBound = storm::utility::one<ValueType>();
    }

    return storm::Interval(lowerBound, upperBound);
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<storm::Interval>> learnIMDPFromMDPByClopperPearsonUntilMaxSamples(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, double lambda, std::size_t numberOfSamples) {
    if constexpr (!storm::IsIntervalType<ValueType>) {
        if (!model->isOfType(storm::models::ModelType::Dtmc) && !model->isOfType(storm::models::ModelType::Mdp)) {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "We only support learning an interval model for DTMCs and MDPs as SUL.");
        }
        STORM_LOG_THROW(lambda > 0.0 && lambda < 1.0, storm::exceptions::InvalidOperationException,
                        "The local successor-distribution failure probability must be in the open interval (0, 1).");
        STORM_PRINT_AND_LOG("Learning Clopper-Pearson intervals with local successor-distribution failure probability "
                            << lambda << " allocated over outgoing successors.\n");

        using IntervalType = storm::Interval;
        using RewardModelType = storm::models::sparse::StandardRewardModel<IntervalType>;

        std::unordered_map<std::size_t, std::uint_fast64_t> choiceSampleCount;
        std::unordered_map<std::size_t, std::unordered_map<storm::storage::sparse::state_type, std::uint_fast64_t>> choiceSuccessorSampleCount;
        storm::storage::BitVector choiceHasSingleSuccessor(model->getNumberOfChoices(), false);

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

                if (successorDistribution.size() == 1) {
                    choiceHasSingleSuccessor.set(currentChoice, true);
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

        // Construct IMDP from the results of the sampling.
        auto const& oldMatrix = model->getTransitionMatrix();

        bool const isDeterministicModel = model->isOfType(storm::models::ModelType::Dtmc);
        storm::storage::SparseMatrixBuilder<IntervalType> builder(oldMatrix.getRowCount(), oldMatrix.getColumnCount(), oldMatrix.getNonzeroEntryCount(), true,
                                                                  !isDeterministicModel, isDeterministicModel ? 0 : oldMatrix.getRowGroupCount());
        ValueType maximumFeasibleL1Diameter = storm::utility::zero<ValueType>();

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            if (!isDeterministicModel) {
                builder.newRowGroup(stateChoiceIndices[s]);
            }

            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);
                auto componentFailureProbability = storm::api::getComponentFailureProbability(lambda, storm::api::countSuccessors(begin, end));
                ValueType lowerSum = storm::utility::zero<ValueType>();
                ValueType upperSum = storm::utility::zero<ValueType>();

                auto successorState = begin;
                while (successorState != end) {
                    IntervalType estimatedInterval;

                    // Small support optimization: if we know that the choice has only one possible successor, we can just set the interval [1, 1].
                    if (choiceHasSingleSuccessor.get(currentChoice)) {
                        estimatedInterval = IntervalType(storm::utility::one<ValueType>(), storm::utility::one<ValueType>());
                    } else {  // Otherwise, we compute the Clopper-Pearson intervals.
                        estimatedInterval =
                            storm::api::getClopperPearsonInterval<ValueType>(choiceSuccessorSampleCount[currentChoice][successorState->getColumn()],
                                                                             choiceSampleCount[currentChoice], componentFailureProbability);
                    }

                    lowerSum += estimatedInterval.lower();
                    upperSum += estimatedInterval.upper();
                    builder.addNextValue(currentChoice, successorState->getColumn(), estimatedInterval);

                    successorState++;

                    STORM_LOG_THROW((successorState == end) || !choiceHasSingleSuccessor.get(currentChoice), storm::exceptions::InvalidOperationException,
                                    "Choice was supposed to only have one successor. How?");
                }
                ValueType lowerSlack = storm::utility::one<ValueType>() - lowerSum;
                ValueType upperSlack = upperSum - storm::utility::one<ValueType>();
                ValueType feasibleMass = std::min(lowerSlack, upperSlack);
                if (feasibleMass < storm::utility::zero<ValueType>()) {
                    feasibleMass = storm::utility::zero<ValueType>();
                }
                ValueType feasibleL1Diameter = static_cast<ValueType>(2.0) * feasibleMass;
                maximumFeasibleL1Diameter = std::max(maximumFeasibleL1Diameter, feasibleL1Diameter);
            }
        }

        STORM_PRINT_AND_LOG("Maximum feasible L1-diameter over all learned choices: " << maximumFeasibleL1Diameter << ".\n");

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

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<storm::Interval>> learnIMDPFromMDPByClopperPearsonUntilL1Width(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, double lambda, double maxL1Width, std::size_t maxNumberOfSamples) {
    if constexpr (!storm::IsIntervalType<ValueType>) {
        if (!model->isOfType(storm::models::ModelType::Dtmc) && !model->isOfType(storm::models::ModelType::Mdp)) {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "We only support learning an interval model for DTMCs and MDPs as SUL.");
        }

        STORM_LOG_THROW(maxL1Width > 0, storm::exceptions::InvalidOperationException, "Maximum L1 Clopper-Pearson interval width must be positive.");
        STORM_LOG_THROW(maxNumberOfSamples > 0, storm::exceptions::InvalidOperationException, "Maximum number of samples must be positive.");
        STORM_LOG_THROW(lambda > 0.0 && lambda < 1.0, storm::exceptions::InvalidOperationException,
                        "The local successor-distribution failure probability must be in the open interval (0, 1).");
        STORM_PRINT_AND_LOG("Learning Clopper-Pearson intervals with local successor-distribution failure probability "
                            << lambda << " allocated over outgoing successors.\n");

        using IntervalType = storm::Interval;
        using RewardModelType = storm::models::sparse::StandardRewardModel<IntervalType>;

        std::unordered_map<std::size_t, std::uint_fast64_t> choiceSampleCount;
        std::unordered_map<std::size_t, std::unordered_map<storm::storage::sparse::state_type, std::uint_fast64_t>> choiceSuccessorSampleCount;
        storm::storage::BitVector choiceHasSingleSuccessor(model->getNumberOfChoices(), false);
        std::uint_fast64_t maximumChoiceSampleCount = 0;

        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

        auto stateChoiceIndices = model->getTransitionMatrix().getRowGroupIndices();

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);

                storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> successorDistribution;
                auto successorState = begin;
                while (successorState != end) {
                    successorDistribution.addProbability(successorState->getColumn(), successorState->getValue());
                    successorState++;
                }

                if (successorDistribution.size() == 1) {
                    choiceHasSingleSuccessor.set(currentChoice, true);
                    continue;
                }
                auto componentFailureProbability = storm::api::getComponentFailureProbability(lambda, successorDistribution.size());

                bool l1WidthTargetReached = false;
                std::uint_fast64_t batchSize = 100;
                while (!l1WidthTargetReached && choiceSampleCount[currentChoice] < maxNumberOfSamples) {
                    auto remainingSamples = maxNumberOfSamples - choiceSampleCount[currentChoice];
                    auto currentBatchSize = std::min<std::uint_fast64_t>(batchSize, remainingSamples);
                    for (std::uint_fast64_t sample = 0; sample < currentBatchSize; ++sample) {
                        auto sampledSuccessorState = successorDistribution.sampleFromDistribution(uniformDistribution(rng));
                        choiceSuccessorSampleCount[currentChoice][sampledSuccessorState] = choiceSuccessorSampleCount[currentChoice][sampledSuccessorState] + 1;
                        choiceSampleCount[currentChoice] = choiceSampleCount[currentChoice] + 1;
                    }

                    ValueType lowerSum = storm::utility::zero<ValueType>();
                    ValueType upperSum = storm::utility::zero<ValueType>();
                    successorState = begin;
                    while (successorState != end) {
                        auto interval = storm::api::getClopperPearsonInterval<ValueType>(choiceSuccessorSampleCount[currentChoice][successorState->getColumn()],
                                                                                         choiceSampleCount[currentChoice], componentFailureProbability);
                        lowerSum += interval.lower();
                        upperSum += interval.upper();
                        successorState++;
                    }
                    ValueType lowerSlack = storm::utility::one<ValueType>() - lowerSum;
                    ValueType upperSlack = upperSum - storm::utility::one<ValueType>();
                    ValueType feasibleMass = std::min(lowerSlack, upperSlack);
                    if (feasibleMass < storm::utility::zero<ValueType>()) {
                        feasibleMass = storm::utility::zero<ValueType>();
                    }
                    ValueType feasibleL1Diameter = static_cast<ValueType>(2.0) * feasibleMass;
                    l1WidthTargetReached = feasibleL1Diameter <= maxL1Width;

                    batchSize = std::min<std::uint_fast64_t>(batchSize * 2, 10000);
                }

                maximumChoiceSampleCount = std::max(maximumChoiceSampleCount, choiceSampleCount[currentChoice]);

                STORM_LOG_THROW(l1WidthTargetReached, storm::exceptions::InvalidOperationException,
                                "Could not learn all state-action Clopper-Pearson feasible L1 diameters below "
                                    << maxL1Width << " within " << maxNumberOfSamples << " samples for choice " << currentChoice
                                    << ". Maximum samples needed so far: " << maximumChoiceSampleCount << ".");
            }
        }

        STORM_PRINT_AND_LOG("Maximum samples needed for L1-width learning over all nontrivial choices: " << maximumChoiceSampleCount << ".\n");

        auto const& oldMatrix = model->getTransitionMatrix();

        bool const isDeterministicModel = model->isOfType(storm::models::ModelType::Dtmc);
        storm::storage::SparseMatrixBuilder<IntervalType> builder(oldMatrix.getRowCount(), oldMatrix.getColumnCount(), oldMatrix.getNonzeroEntryCount(), true,
                                                                  !isDeterministicModel, isDeterministicModel ? 0 : oldMatrix.getRowGroupCount());

        for (std::size_t s = 0; s < model->getNumberOfStates(); ++s) {
            if (!isDeterministicModel) {
                builder.newRowGroup(stateChoiceIndices[s]);
            }

            auto numberOfChoices = stateChoiceIndices[s + 1] - stateChoiceIndices[s];
            for (auto choiceOffset = 0; choiceOffset < numberOfChoices; choiceOffset++) {
                auto currentChoice = stateChoiceIndices[s] + choiceOffset;
                auto begin = model->getTransitionMatrix().begin(currentChoice), end = model->getTransitionMatrix().end(currentChoice);
                auto componentFailureProbability = storm::api::getComponentFailureProbability(lambda, storm::api::countSuccessors(begin, end));

                auto successorState = begin;
                while (successorState != end) {
                    IntervalType estimatedInterval;

                    if (choiceHasSingleSuccessor.get(currentChoice)) {
                        estimatedInterval = IntervalType(storm::utility::one<ValueType>(), storm::utility::one<ValueType>());
                    } else {
                        estimatedInterval =
                            storm::api::getClopperPearsonInterval<ValueType>(choiceSuccessorSampleCount[currentChoice][successorState->getColumn()],
                                                                             choiceSampleCount[currentChoice], componentFailureProbability);
                    }

                    builder.addNextValue(currentChoice, successorState->getColumn(), estimatedInterval);

                    successorState++;

                    STORM_LOG_THROW((successorState == end) || !choiceHasSingleSuccessor.get(currentChoice), storm::exceptions::InvalidOperationException,
                                    "Choice was supposed to only have one successor. How?");
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

        if (isDeterministicModel) {
            return std::make_shared<storm::models::sparse::Dtmc<IntervalType>>(std::move(components));
        }

        return std::make_shared<storm::models::sparse::Mdp<IntervalType>>(std::move(components));
    }
}

}  // namespace api
}  // namespace storm

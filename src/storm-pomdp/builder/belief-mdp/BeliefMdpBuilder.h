#pragma once

#include "storm-pomdp/builder/belief-mdp/BeliefExploration.h"
#include "storm-pomdp/builder/belief-mdp/BeliefExplorationHeuristic.h"
#include "storm/logic/Formulas.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/sparse/ModelComponents.h"

#include "storm/exceptions/UnexpectedException.h"

namespace storm::pomdp::builder {

struct BeliefMdpPropertyInformation {
    enum class Kind { ReachabilityProbability, ExpectedTotalReachabilityReward };
    Kind kind;
    std::set<uint32_t> targetObservations;
    std::optional<std::string> rewardModelName;
    storm::OptimizationDirection dir;
};

template<BeliefExplorationMode Mode, typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
class BeliefMdpBuilder {
   public:
    using BeliefId = storm::pomdp::beliefs::BeliefId;
    using StateId = storm::pomdp::beliefs::BeliefStateType;
    using PomdpValueType = typename PomdpType::ValueType;
    using BeliefValueType = typename BeliefType::ValueType;
    using ObservationType = uint32_t;
    using ExplorationInformation = typename BeliefExploration<Mode, BeliefMdpValueType, PomdpType, BeliefType>::ExplorationInformation;

    BeliefMdpBuilder(PomdpType const& pomdp, BeliefMdpPropertyInformation propertyInformation) : explorer(pomdp), propertyInformation(propertyInformation) {
        if (propertyInformation.rewardModelName) {
            explorer.setRewardModelForObjective(*propertyInformation.rewardModelName);
        }
    }

    ExplorationInformation explore(/* explorationHeuristic, termination, abstraction*/) {
        BeliefExplorationHeuristic h;
        auto info = explorer.initializeExploration(h);
        explorer.performExploration(
            info, h, []() { return false; }, [this](BeliefType const& bel) { return propertyInformation.targetObservations.count(bel.observation()) != 0; },
            storm::pomdp::beliefs::NoAbstraction());
    }

    std::shared_ptr<storm::models::sparse::Mdp<BeliefMdpValueType>> build(ExplorationInformation const& info, /* CutOffData */) {
        // TODO: reachability analysis?
        return createMdpFromExplorationInfo(info);
    }

    std::shared_ptr<storm::logic::Formula const> createFormulaForBeliefMdp() {
        STORM_LOG_ASSERT(propertyInformation.kind == BeliefMdpPropertyInformation::Kind::ReachabilityProbability ||
                             propertyInformation.kind == BeliefMdpPropertyInformation::Kind::ExpectedTotalReachabilityReward,
                         "Unexpected kind of property.");
        switch (propertyInformation.kind) {
            case BeliefMdpPropertyInformation::Kind::ReachabilityProbability: {
                auto target = std::make_shared<storm::logic::AtomicLabelFormula const>("target");
                auto eventuallyTarget = std::make_shared<storm::logic::EventuallyFormula const>(target, storm::logic::FormulaContext::Probability);
                return std::make_shared<storm::logic::ProbabilityOperatorFormula const>(eventuallyTarget,
                                                                                        storm::logic::OperatorInformation(propertyInformation.dir));
            }
            case BeliefMdpPropertyInformation::Kind::ExpectedTotalReachabilityReward: {
                auto bottom = std::make_shared<storm::logic::AtomicLabelFormula const>("bottom");
                auto eventuallyBottom = std::make_shared<storm::logic::EventuallyFormula const>(bottom, storm::logic::FormulaContext::Reward,
                                                                                                storm::logic::RewardAccumulation(true, false, false));
                return std::make_shared<storm::logic::RewardOperatorFormula const>(eventuallyBottom, propertyInformation.rewardModelName.value(),
                                                                                   storm::logic::OperatorInformation(propertyInformation.dir));
            }
        }
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unhandled case.");
    }

   private:
    std::shared_ptr<storm::models::sparse::Mdp<BeliefMdpValueType>> createMdpFromExplorationInfo(ExplorationInformation const& info) {
        STORM_LOG_ASSERT(propertyInformation.kind == BeliefMdpPropertyInformation::Kind::ReachabilityProbability ||
                             propertyInformation.kind == BeliefMdpPropertyInformation::Kind::ExpectedTotalReachabilityReward,
                         "Unexpected kind of property.");
        bool const reachabilityProbability = propertyInformation.kind == BeliefMdpPropertyInformation::Kind::ReachabilityProbability;
        uint64_t const numExtraStates = reachabilityProbability ? 2ull : 1ull;
        uint64_t const numStates = info.matrix.groups() + numExtraStates;
        uint64_t const numChoices = info.matrix.rows() + numExtraStates;
        uint64_t const targetState = numStates - numExtraStates;
        uint64_t const bottomState = numStates - 1;
        std::vector<BeliefMdpValueType> actionRewards;
        if (!reachabilityProbability) {
            actionRewards.reserve(numChoices);
            actionRewards.insert(actionRewards.end(), info.actionRewards.begin(), info.actionRewards.end());
            actionRewards.push_back(storm::utility::zero<BeliefMdpValueType>());
            STORM_LOG_ASSERT(numChoices == actionRewards.size(),
                             "Unexpected size of action rewards: Expected " << numChoices << " got " << actionRewards.size() << ".");
        }
        storm::storage::SparseMatrixBuilder<BeliefMdpValueType> transitionBuilder(numChoices, numStates, 0, true, true, numStates);
        for (uint64_t state = 0; state < numStates - numExtraStates; ++state) {
            uint64_t choice = info.matrix.rowGroupIndices[state];
            transitionBuilder.newRowGroup(choice);
            for (uint64_t const groupEnd = info.matrix.rowGroupIndices[state + 1]; choice < groupEnd; ++choice) {
                BeliefMdpValueType probabilityToBottom = storm::utility::zero<BeliefMdpValueType>();
                BeliefMdpValueType probabilityToTarget = storm::utility::zero<BeliefMdpValueType>();
                for (uint64_t entryIndex = info.matrix.rowIndications[choice]; entryIndex < info.matrix.rowIndications[choice + 1]; ++entryIndex) {
                    auto const& entry = info.matrix.transitions[entryIndex];
                    if (auto findIt = info.exploredBeliefs.find(entry.targetBelief); findIt != info.exploredBeliefs.end()) {
                        // Transition to explored belief
                        transitionBuilder.addNextValue(choice, findIt->second, entry.probability);
                    } else {
                        // Transition to unexplored belief
                        auto const& successorBelief = info.collectedBeliefs.getBeliefFromId(entry.targetBelief);
                        bool const successorIsTarget = propertyInformation.targetObservations.count(successorBelief.observation()) != 0;
                        if (reachabilityProbability) {
                            if (successorIsTarget) {
                                probabilityToTarget += entry.probability;
                            } else {
                                auto const cutOffValue = entry.probability * computeCutOffValue(successorBelief);
                                probabilityToTarget += cutOffValue;
                                probabilityToBottom += entry.probability - cutOffValue;
                            }
                        } else {
                            probabilityToBottom += entry.probability;
                            if (!successorIsTarget) {
                                actionRewards[choice] += entry.probability * computeCutOffValue(successorBelief);
                            }
                        }
                    }
                }
                if (reachabilityProbability && !storm::utility::isZero(probabilityToTarget)) {
                    transitionBuilder.addNextValue(choice, targetState, probabilityToTarget);
                }
                if (!storm::utility::isZero(probabilityToBottom)) {
                    transitionBuilder.addNextValue(choice, bottomState, probabilityToBottom);
                }
            }
        }
        if (reachabilityProbability) {
            transitionBuilder.newRowGroup(numChoices - 2);
            transitionBuilder.addNextValue(numChoices - 2, targetState, storm::utility::one<BeliefMdpValueType>());
        }
        transitionBuilder.newRowGroup(numChoices - 1);
        transitionBuilder.addNextValue(numChoices - 1, bottomState, storm::utility::one<BeliefMdpValueType>());

        storm::models::sparse::StateLabeling stateLabeling(numStates);
        stateLabeling.addLabel("bottom");
        stateLabeling.addLabelToState("bottom", bottomState);
        stateLabeling.addLabel("init");
        stateLabeling.addLabelToState("init", info.exploredBeliefs.at(info.initialBelief));

        if (reachabilityProbability) {
            stateLabeling.addLabel("target");
            stateLabeling.addLabelToState("target", targetState);
        }
        storm::storage::sparse::ModelComponents<BeliefMdpValueType> components(transitionBuilder.build(), std::move(stateLabeling));

        if (!reachabilityProbability) {
            storm::models::sparse::StandardRewardModel<BeliefMdpValueType> rewardModel(std::nullopt, std::move(actionRewards));
            components.rewardModels.emplace(propertyInformation.rewardModelName.value(), std::move(rewardModel));
        }

        return std::make_shared<storm::models::sparse::Mdp<BeliefMdpValueType>>(std::move(components));
    }

    BeliefMdpValueType computeCutOffValue(BeliefType const& belief) {
        // TODO
        assert(false);
    }
    BeliefExploration<Mode, BeliefMdpValueType, PomdpType, BeliefType> explorer;
    BeliefMdpPropertyInformation propertyInformation;
};
}  // namespace storm::pomdp::builder

#include "SparseModelDARewardProduct.h"

#include <storm/environment/modelchecker/ModelCheckerEnvironment.h>

#include <boost/spirit/home/qi/string.hpp>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/environment/Environment.h"
#include "storm/transformer/DARewardProductBuilder.h"

namespace storm {
namespace storage {

template<typename ValueType, typename RewardModelType>
void printMDP(
    const storm::storage::SparseMatrix<ValueType>& transitionMatrix,
    const storm::models::sparse::StateLabeling& stateLabeling,
    const std::unordered_map<std::string, RewardModelType>& rewardModels,
    int numStatesToPrint=0,
    bool printCompleteRow=true
) {
    std::cout << "Markov Decision Process (MDP) Details:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    if (!numStatesToPrint) numStatesToPrint = transitionMatrix.getRowGroupCount();

    // 1. Zustandsbeschreibung (Labels)
    std::cout << "State Labels:" << std::endl;
    for (uint64_t state = 0; state < numStatesToPrint; ++state) {
        std::cout << "  State " << state << ": ";
        auto labels = stateLabeling.getLabelsOfState(state);
        if (!labels.empty()) {
            for (const auto& label : labels) {
                std::cout << label << " ";
            }
        } else {
            std::cout << "(No labels)";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // 2. Transitionsmatrix (Wahrscheinlichkeiten)
    std::cout << "Transition Matrix:" << std::endl;
    for (uint64_t state = 0; state < numStatesToPrint; ++state) {
        std::cout << "  State " << state << ":" << std::endl;

        for (auto const& choice : transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);

            std::cout << "    Choice " << choice << ": ";
            uint64_t currentColumn = 0;

            if (printCompleteRow) {
                for (auto const& entry : row) {
                    for (; currentColumn < entry.getColumn(); ++currentColumn) {
                        std::cout << "0 ";
                    }
                    std::cout << entry.getValue() << " ";
                    ++currentColumn;
                }

                for (; currentColumn < numStatesToPrint; ++currentColumn) {
                    std::cout << "0 ";
                }
            } else {
                for (auto const& entry : row) {
                    std::cout << "(" << entry.getColumn() << ": " << entry.getValue() << ") ";
                }
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    // 3. Belohnungsmodelle
    std::cout << "Reward Models:" << std::endl;
    for (const auto& [rewardName, rewardModel] : rewardModels) {
        std::cout << "  Reward Model \"" << rewardName << "\":" << std::endl;
        if (rewardModel.hasStateRewards()) {
            for (uint64_t state = 0; state < numStatesToPrint; ++state) {
                auto reward = rewardModel.getStateReward(state);
                std::cout << "    State " << state << ": " << reward << std::endl;
            }
        }
        if (rewardModel.hasStateActionRewards()) {
            uint64_t lastChoice = transitionMatrix.getRowGroupIndices(numStatesToPrint-1).back();
            for (uint64_t choice = 0; choice < lastChoice; ++choice) {
                auto reward = rewardModel.getStateActionReward(choice);
                std::cout << "    Choice " << choice << ": " << reward << std::endl;
            }
        }
    }

    std::cout << "-------------------------------------" << std::endl;
}


template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Mdp<ValueType, RewardModelType>> SparseModelDARewardProduct<ValueType, RewardModelType>::build() {
    storm::storage::BitVector initialStatesProduct;
    //printMDP(originalModel.getTransitionMatrix(), originalModel.getStateLabeling(), originalModel.getRewardModels());
    auto [productModel, acceptanceConditions, indexToModelState] = modelchecker::helper::SparseLTLHelper<ValueType, true>::buildFromFormulas(originalModel, formulas, env);
    //printMDP(productModel.getTransitionMatrix(), productModel.getStateLabeling(), productModel.getRewardModels());

    if (env.modelchecker().isLtl2daToolSet()) {
        auto rewardModels = buildRewardModelsForLDBA(productModel, acceptanceConditions, indexToModelState);
        //printMDP(productModel.getTransitionMatrix(), productModel.getStateLabeling(), rewardModels);
        return std::make_shared<Mdp>(productModel.getTransitionMatrix(), productModel.getStateLabeling(), rewardModels);
    }

    transformer::DARewardProductBuilder<ValueType, RewardModelType> builder(productModel, acceptanceConditions, indexToModelState, originalModel);
    auto result = builder.build();

    auto rewardModels = buildRewardModels(result->getTransitionMatrix(), result->getStateToModelState(), result->getActionToModelAction());
    auto LTLRewardModel = buildLTLRewardModel(result->getTransitionMatrix(), result->getReachingAccEcChoices());
    rewardModels.merge(LTLRewardModel);

    auto stateLabeling = buildStateLabeling(result->getTransitionMatrix(), result->getStateToModelState(), result->getInitialStates());

    //printMDP(result->getTransitionMatrix(), stateLabeling, rewardModels);

    return std::make_shared<Mdp>(result->getTransitionMatrix(), stateLabeling, rewardModels);
}

template<typename ValueType, typename RewardModelType>
std::unordered_map<std::string, RewardModelType> SparseModelDARewardProduct<ValueType, RewardModelType>::buildLTLRewardModel(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<std::list<uint64_t>> const& reachingAccECsChoices) {
    typedef typename RewardModelType::ValueType RewardValueType;
    std::unordered_map<std::string, RewardModelType> result;
    const uint64_t numberOfLTLObjectives = std::log2(reachingAccECsChoices.size());

    for (int i = 0; i < numberOfLTLObjectives; i++) {
        // add rewards for reaching an accepting end component
        std::vector<RewardValueType> stateActionRewards(resultTransitionMatrix.getRowCount(), storm::utility::zero<RewardValueType>());
        for (int j = 0; j < reachingAccECsChoices.size(); j++) {
            if (!(j & (1 << i))) continue;

            for (auto const& choice: reachingAccECsChoices[j]) {
                stateActionRewards[choice] = 1;
            }
        }

        std::string name = "accEc_" + std::to_string(i);
        result.insert(std::make_pair(name, RewardModelType(std::nullopt, stateActionRewards)));
    }

    return result;
}

template<typename ValueType, typename RewardModelType>
std::unordered_map<std::string, RewardModelType> SparseModelDARewardProduct<ValueType, RewardModelType>::buildRewardModels(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, std::vector<uint64_t> const& choiceToModelChoice) {
    typedef typename RewardModelType::ValueType RewardValueType;
    std::unordered_map<std::string, RewardModelType> result;
    uint64_t numResStates = resultTransitionMatrix.getRowGroupCount();
    uint64_t numResChoices = resultTransitionMatrix.getRowCount();

    for (auto const rewardModel: originalModel.getRewardModels()) {
        std::optional<std::vector<RewardValueType>> stateRewards;
        if (rewardModel.second.hasStateRewards()) {
            stateRewards = std::vector<RewardValueType>(numResStates, storm::utility::zero<RewardValueType>());
            for (uint64_t state = 0; state < numResStates; ++state) {
                if (stateToModelState[state] == std::numeric_limits<uint64_t>::max()) continue;

                stateRewards.value()[state] = rewardModel.second.getStateReward(stateToModelState[state]);
            }
        }

        std::optional<std::vector<RewardValueType>> actionRewards;
        if (rewardModel.second.hasStateActionRewards()) {
            actionRewards = std::vector<RewardValueType>(numResChoices, storm::utility::zero<RewardValueType>());

            for (uint64_t choice = 0; choice < numResChoices; ++choice) {
                if (choiceToModelChoice[choice] == std::numeric_limits<uint64_t>::max()) continue;

                actionRewards.value()[choice] += rewardModel.second.getStateActionReward(choiceToModelChoice[choice]);
            }
        }

        result.insert(std::make_pair(rewardModel.first, RewardModelType(stateRewards, actionRewards)));
    }

    return result;
}

template<typename ValueType, typename RewardModelType>
storm::models::sparse::StateLabeling SparseModelDARewardProduct<ValueType, RewardModelType>::buildStateLabeling(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, storm::storage::BitVector const& initialStates) {
    uint64_t modelStateCount = originalModel.getNumberOfStates();
    uint64_t numResStates = resultTransitionMatrix.getRowGroupCount();
    storm::models::sparse::StateLabeling resultLabeling(numResStates);

    for (auto const& label: originalModel.getStateLabeling().getLabels()) {
        resultLabeling.addLabel(label);
    }

    for (uint64_t state = 0; state < numResStates; ++state) {
        uint64_t modelState = stateToModelState[state];
        if (modelState == std::numeric_limits<uint64_t>::max()) {
            continue;
        }

        for (auto const& label : originalModel.getLabelsOfState(modelState)) {
            if (label != "init") {
                resultLabeling.addLabelToState(label, state);
            }
        }
    }

    resultLabeling.setStates("init", initialStates);

    return resultLabeling;
}

template<typename ValueType, typename RewardModelType>
std::unordered_map<std::string, RewardModelType> SparseModelDARewardProduct<ValueType, RewardModelType>::buildRewardModelsForLDBA(
    Mdp& productModel, std::vector<storm::automata::AcceptanceCondition::ptr> const& acceptanceConditions, std::vector<uint64_t> indexToModelState) {
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(productModel.getTransitionMatrix(), productModel.getBackwardTransitions());
    typedef typename RewardModelType::ValueType RewardValueType;
    std::unordered_map<std::string, RewardModelType> rewardModels;
    uint64_t numStates = productModel.getTransitionMatrix().getRowGroupCount();
    uint64_t numChoices = productModel.getTransitionMatrix().getRowCount();
    storm::storage::BitVector acc_states(numStates);
    storm::storage::BitVector initial_states(numStates);
    for (auto const& state: mecs[0].getStateSet()) {
        initial_states.set(state);
    }

    for (int i = 0; i < acceptanceConditions.size(); i++) {
        std::optional<std::vector<RewardValueType>> stateRewards = std::vector<RewardValueType>(numStates, storm::utility::zero<RewardValueType>());
        auto expr = acceptanceConditions[i]->getAcceptanceExpression();
        STORM_LOG_ASSERT(expr->isAtom() && expr->getAtom().getType() == cpphoafparser::AtomAcceptance::TEMPORAL_INF,
                         "For Büchi acceptance, the acceptance expression has to be of the form INF(0)");

        auto const& acceptanceSet = expr->getAtom().getAcceptanceSet();
        auto const& acceptingStates = acceptanceConditions[i]->getAcceptanceSet(acceptanceSet);

        for (auto const& ec : mecs) {
            //if (ec.size() != 426) continue;
            bool containsAcceptingStates = false;
            for (auto const& state : acceptingStates) {
                if (ec.containsState(state)) {
                    containsAcceptingStates = true;
                    break;
                }
            }
            if (!containsAcceptingStates) continue;

            for (auto const& state : ec.getStateSet()) {
                stateRewards.value()[state] = 1;
                acc_states.set(state);
            }
        }
        rewardModels.insert(std::make_pair("accEc_" + std::to_string(i), RewardModelType(stateRewards)));
    }

    for (auto const rewardModel : originalModel.getRewardModels()) {
        std::optional<std::vector<RewardValueType>> stateRewards;
        if (rewardModel.second.hasStateRewards()) {
            stateRewards = std::vector<RewardValueType>(numStates, storm::utility::zero<RewardValueType>());
            for (uint64_t state = 0; state < numStates; ++state) {
                //if (!initial_states[state]) continue; //!acc_states[state] &&
                stateRewards.value()[state] = rewardModel.second.getStateReward(indexToModelState[state]);
            }
        }

        /* todo later
        std::optional<std::vector<RewardValueType>> actionRewards;
        if (rewardModel.second.hasStateActionRewards()) {
            actionRewards = std::vector<RewardValueType>(numChoices, storm::utility::zero<RewardValueType>());
            for (uint64_t state = 0; state < originalModel.getNumberOfStates(); ++state) {
                for (uint64_t choice = 0; choice < originalModel.getTransitionMatrix().getRowGroupSize(state); ++choice) {
                    auto modelChoice = BitVector(originalModel.getTransitionMatrix().getRowCount(), false);
                    modelChoice.set(originalModel.getTransitionMatrix().getRowGroupIndices(state)[choice], true);
                    auto liftedChoices = productModel.liftFromModel(modelChoice);

                    for (auto const& productChoice : liftedChoices) {
                        actionRewards.value()[productChoice] =
                            rewardModel.second.getStateActionReward(originalModel.getTransitionMatrix().getRowGroupIndices(state)[choice]);
                    }
                }
            }
        }
        */
        rewardModels.insert(std::make_pair(rewardModel.first, RewardModelType(stateRewards)));
    }

    return rewardModels;
}

template class SparseModelDARewardProduct<double, storm::models::sparse::StandardRewardModel<double>>;

template class SparseModelDARewardProduct<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}
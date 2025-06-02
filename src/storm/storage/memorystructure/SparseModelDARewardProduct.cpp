
#include "SparseModelDARewardProduct.h"

#include <boost/spirit/home/qi/string.hpp>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/transformer/DARewardProductBuilder.h"

namespace storm {
namespace storage {

template<typename ValueType, typename RewardModelType>
void printMDP(
    const storm::storage::SparseMatrix<ValueType>& transitionMatrix,
    const storm::models::sparse::StateLabeling& stateLabeling,
    const std::unordered_map<std::string, RewardModelType>& rewardModels
) {
    std::cout << "Markov Decision Process (MDP) Details:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // 1. Zustandsbeschreibung (Labels)
    std::cout << "State Labels:" << std::endl;
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
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
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
        std::cout << "  State " << state << ":" << std::endl;

        for (auto const& choice : transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);

            std::cout << "    Choice " << choice << ": ";
            uint64_t currentColumn = 0;

            for (auto const& entry : row) {
                // Fülle Lücken mit Nullen bis zur nächsten belegten Spalte
                for (; currentColumn < entry.getColumn(); ++currentColumn) {
                    std::cout << "0 ";
                }
                // Gebe Wert aus
                std::cout << entry.getValue() << " ";
                ++currentColumn;
            }

            // Fülle verbleibende Nullspalten auf
            for (; currentColumn < transitionMatrix.getColumnCount(); ++currentColumn) {
                std::cout << "0 ";
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
            for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
                auto reward = rewardModel.getStateReward(state);
                std::cout << "    State " << state << ": " << reward << std::endl;
            }
        }
        if (rewardModel.hasStateActionRewards()) {
            for (uint64_t choice = 0; choice < transitionMatrix.getRowCount(); ++choice) {
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
    //std::tie(product, initialStatesProduct)
    auto [productModel, acceptanceConditions, indexToModelState] = modelchecker::helper::SparseLTLHelper<ValueType, true>::buildFromFormulas(originalModel, formulas);

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

    for (auto const rewardModel: originalModel.getRewardModels()) {
        std::optional<std::vector<RewardValueType>> stateRewards;
        if (rewardModel.second.hasStateRewards()) {
            stateRewards = std::vector<RewardValueType>(numResStates, storm::utility::zero<RewardValueType>());
            for (uint64_t state = 0; state < numResStates; ++state) {
                if (stateToModelState[state] != std::numeric_limits<uint64_t>::max()) {
                    stateRewards.value()[state] = rewardModel.second.getStateReward(stateToModelState[state]);
                }
            }
        }

        std::optional<std::vector<RewardValueType>> stateActionRewards;
        if (rewardModel.second.hasStateActionRewards()) {
            stateActionRewards = std::vector<RewardValueType>(numResStates, storm::utility::zero<RewardValueType>());

            for (uint64_t choice = 0; choice < resultTransitionMatrix.getRowCount(); ++choice) {
                if (choiceToModelChoice[choice] != std::numeric_limits<uint64_t>::max()) {
                    stateActionRewards.value()[choice] += rewardModel.second.getStateActionReward(choiceToModelChoice[choice]);
                }
            }
        }

        result.insert(std::make_pair(rewardModel.first, RewardModelType(stateRewards, stateActionRewards)));
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

    //std::cout << initialStates << std::endl;
    resultLabeling.setStates("init", initialStates);

    return resultLabeling;
}

template class SparseModelDARewardProduct<double, storm::models::sparse::StandardRewardModel<double>>;

template class SparseModelDARewardProduct<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}
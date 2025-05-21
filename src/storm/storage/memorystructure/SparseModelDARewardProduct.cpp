
#include "SparseModelDARewardProduct.h"
#include "storm/transformer/DARewardProductBuilder.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm {
namespace storage {

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Mdp<ValueType, RewardModelType>> SparseModelDARewardProduct<ValueType, RewardModelType>::build() {
    storm::storage::BitVector initialStatesProduct;
    //std::tie(product, initialStatesProduct)
    typename transformer::DAMultiProduct<Mdp>::ptr test = modelchecker::helper::SparseLTLHelper<ValueType, true>::buildFromFormulas(originalModel, formulas);

    /*transformer::DARewardProductBuilder<ValueType, RewardModelType> builder(*product, originalModel, initialStatesProduct);
    auto result = builder.build();

    auto rewardModels = buildRewardModels(result->getTransitionMatrix(), result->getStateToModelState(), result->getActionToModelAction(), result->getReachingAccEcChoices());
    auto stateLabeling = buildStateLabeling(result->getTransitionMatrix(), result->getStateToModelState(), result->getInitialStates());

    return std::make_shared<Mdp>(result->getTransitionMatrix(), stateLabeling, rewardModels);
    */
    return nullptr;
}

template<typename ValueType, typename RewardModelType>
std::unordered_map<std::string, RewardModelType> SparseModelDARewardProduct<ValueType, RewardModelType>::buildRewardModels(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, std::vector<uint64_t> const& choiceToModelChoice, std::list<uint64_t> const& reachingAccECsChoices) {
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

    // add rewards for reaching an accepting end component
    std::vector<RewardValueType> stateActionRewards(resultTransitionMatrix.getRowCount(), storm::utility::zero<RewardValueType>());
    for (auto const& choice: reachingAccECsChoices) {
        stateActionRewards[choice] = 1;
    }
    result.insert(std::make_pair("accEc", RewardModelType(std::nullopt, stateActionRewards)));

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

template class SparseModelDARewardProduct<double, storm::models::sparse::StandardRewardModel<double>>;

template class SparseModelDARewardProduct<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}
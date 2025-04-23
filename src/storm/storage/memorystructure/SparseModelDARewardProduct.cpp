
#include "SparseModelDARewardProduct.h"
#include "storm/transformer/DARewardProductBuilder.h"
#include <gmpxx.h>  //???

namespace storm {
namespace storage {

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Mdp<ValueType, RewardModelType>> SparseModelDARewardProduct<ValueType, RewardModelType>::build() {

    transformer::DARewardProductBuilder<ValueType, RewardModelType> builder(*product);
    auto result = builder.build();
    /*
    auto transitionMatrix = result.transitionMatrix;
    std::cout << "Resulting transition matrix:" << std::endl;
    printMatrix(transitionMatrix);
    STORM_LOG_ASSERT(transitionMatrix.isProbabilistic(), "The resulting transition matrix is not probabilistic");
    std::cout << "Done building the demerged matrix." << std::endl;

    auto rewardModels = buildRewardModels(transitionMatrix, result.stateToModelState, result.reachingAccEcChoices);
    auto stateLabeling = buildStateLabeling(transitionMatrix, stateToModelState, result.initialStates);
    std::cout << "Done building state labeling" << std::endl;

    auto modifiedModel = Mdp(transitionMatrix, stateLabeling, rewardModels);
    return std::make_shared<Mdp>(modifiedModel);
    */
}
/*
template<typename ValueType, typename RewardModelType>
std::unordered_map<std::string, RewardModelType> DARewardProductBuilder<ValueType, RewardModelType>::buildRewardModels(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, std::list<uint64_t> const& reachingAccECsChoices) {
    typedef typename RewardModelType::ValueType RewardValueType;
    std::unordered_map<std::string, RewardModelType> result;
    uint64_t numResStates = resultTransitionMatrix.getRowGroupCount();

    for (auto const rewardModel: model.getRewardModels()) {
        std::optional<std::vector<RewardValueType>> stateRewards;
        if (rewardModel.second.hasStateRewards()) {
            stateRewards = std::vector<RewardValueType>(numResStates, storm::utility::zero<RewardValueType>());
            for (uint64_t state = 0; state < numResStates; ++state) {
                if (stateToModelState[state] != InvalidIndex) {
                    stateRewards.value()[state] = rewardModel.second.getStateReward(stateToModelState[state]);
                }
            }
        }

        result.insert(std::make_pair(rewardModel.first, RewardModelType(stateRewards)));
    }

    // add rewards for reaching an accepting end component
    std::vector<RewardValueType> stateActionRewards(resultTransitionMatrix.getRowGroupCount(), storm::utility::zero<RewardValueType>());
    for (auto const& choice: reachingAccECsChoices) {
        stateActionRewards[choice] = 1;
    }
    result.insert(std::make_pair("accEc", RewardModelType(stateActionRewards)));

    return result;
}

template<typename ValueType, typename RewardModelType>
storm::models::sparse::StateLabeling DARewardProductBuilder<ValueType, RewardModelType>::buildStateLabeling(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, storm::storage::BitVector const& initialStates) {
    uint64_t modelStateCount = model.getNumberOfStates();
    uint64_t numResStates = resultTransitionMatrix.getRowGroupCount();
    storm::models::sparse::StateLabeling resultLabeling(numResStates);

    for (auto const& label: model.getStateLabeling().getLabels()) {
        resultLabeling.addLabel(label);
    }

    for (uint64_t state = 0; state < numResStates; ++state) {
        uint64_t modelState = stateToModelState[state];
        if (modelState == InvalidIndex) {
            continue;
        }

        for (auto const& label : model.getLabelsOfState(modelState)) {
            if (label != "init") {
                resultLabeling.addLabelToState(label, state);
            }
        }
    }
    // brauchen wir memory labels?

    return resultLabeling;
}
*/

template class SparseModelDARewardProduct<double, storm::models::sparse::StandardRewardModel<double>>;

template class SparseModelDARewardProduct<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}
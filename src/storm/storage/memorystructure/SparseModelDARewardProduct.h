#pragma once

#include "storm/logic/Formula.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/transformer/DAProductBuilder.h"
#include "storm/transformer/DARewardProduct.h"
#include "storm/modelchecker/helper/ltl/SparseLTLHelper.h"

namespace storm {
namespace storage {

template <typename ValueType, typename RewardModelType>
class SparseModelDARewardProduct {
public:
    using Mdp = storm::models::sparse::Mdp<ValueType, RewardModelType>;

    SparseModelDARewardProduct(Mdp const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas): originalModel(model), formulas(formulas) {}

    /*!
     * Assembles the modified model by lifting the state labeling and reward models
     * @return a pointer to the modified model
     */
    std::shared_ptr<Mdp> build();

    private:
        Mdp originalModel;
        std::shared_ptr<transformer::DAProduct<Mdp>> product;
        std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas;

        /*!
         *
         * @param transitionMatrix the transition matrix of the modified product model
         * @param stateToProductState mapping of states in the modified to ones in the product model
         * @param initialStates the initial states of the modified model
         * @return the state labeling for the modified model
         */
        storm::models::sparse::StateLabeling buildStateLabeling(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, std::vector<uint64_t> const& stateToProductState, storm::storage::BitVector const& initialStates);

        /*!
         *
         * @param resultTransitionMatrix the transition matrix of the modified model
         * @param stateToProductState mapping of states in the modified to ones in the product model
         * @param reachingAccECsChoices choices that a.s. reach a copy of an accepting ec
         * @return the reward models for the modified model
         */
        std::unordered_map<std::string, RewardModelType> buildRewardModels(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToModelState, std::vector<uint64_t> const& choiceToModelChoice, std::list<uint64_t> const& reachingAccECsChoices);
};
}
}



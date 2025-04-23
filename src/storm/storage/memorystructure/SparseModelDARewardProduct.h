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
    using CheckFormulaCallback = std::function<storm::storage::BitVector(storm::logic::Formula const&)>;
    using Mdp = storm::models::sparse::Mdp<ValueType, RewardModelType>;

    SparseModelDARewardProduct(Mdp const& model, storm::logic::PathFormula const& formula, CheckFormulaCallback const& formulaChecker): originalModel(model), formula(formula), product(modelchecker::helper::SparseLTLHelper<ValueType, true>::buildFromFormula(model, formula, formulaChecker)) {}

    std::shared_ptr<Mdp> build();

    private:
        Mdp originalModel;
        std::shared_ptr<transformer::DAProduct<Mdp>> product;
        storm::logic::PathFormula const& formula;


        /*!
         * Builds the transition matrix of the demerged product MDP.
         * @param acceptance the acceptance condition
         * @param transitionMatrix the transition matrix of the product
         * @param backwardTransitions the reversed transition relation
         * @return A pair consisting of (1) a map from each state in the new model to the state it was copied from and (2) the transition matrix of the demerged model
         */
        storm::transformer::DARewardProduct<ValueType> buildTransitionMatrix(automata::AcceptanceCondition const& acceptance, storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::SparseMatrix<ValueType> const& backwardTransitions);

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
    std::unordered_map<std::string, RewardModelType> buildRewardModels(storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix, std::vector<uint64_t> const& stateToProductState, std::list<uint64_t> const& reachingAccECsChoices);
};
}
}



#pragma once

#include "storm/transformer/DAProductBuilder.h"
#include "storm/logic/Formula.h"
#include "storm/models/sparse/Mdp.h"
//#include "storm/transformer/Product.h"
//#include "storm/transformer/ProductBuilder.h"

namespace storm {
namespace transformer {

template<typename ValueType, typename RewardModelType>
class DARewardProductBuilder {
    public:
        using CheckFormulaCallback = std::function<storm::storage::BitVector(storm::logic::Formula const&)>;
        using Mdp = storm::models::sparse::Mdp<ValueType, RewardModelType>;

        DARewardProductBuilder(Mdp const& model, storm::logic::PathFormula const& formula, CheckFormulaCallback const& formulaChecker): model(model), formula(formula), formulaChecker(formulaChecker) {}

        void build();

        //DARewardProductBuilder(const storm::automata::DeterministicAutomaton &da, const std::vector<storm::storage::BitVector>& statesForAP)
        //    : da(da), statesForAP(statesForAP) {}

        //DAProduct<Mdp<ValueType>::ptr build(const storm::storage::SparseMatrix<ValueType>& originalMatrix, const storm::storage::BitVector& statesOfInterest);

    private:
        storm::logic::PathFormula const& formula;
        Mdp const& model;
        CheckFormulaCallback const& formulaChecker;

        /*!
         * Computes the product MDP for the formula and the model
         * @return the product MDP
         */
        std::shared_ptr<DAProduct<Mdp>> buildProductMDP();

        /*!
         * Computes the product a given determistic automaton and the model
         * @param da deterministic automaton for the LTL formula
         * @param apSatSets the sets where the Aps are satisfied
         * @return the product MDP
         */
        std::shared_ptr<DAProduct<Mdp>> buildDAProduct(storm::automata::DeterministicAutomaton const& da, std::map<std::string, storm::storage::BitVector>& apSatSets);

        /*!
         * Builds the transition matrix of the demerged product MDP.
         * @param acceptance the acceptance condition
         * @param transitionMatrix the transition matrix of the product
         * @param backwardTransitions the reversed transition relation
         * @return A pair consisting of (1) a map from each state in the new model to the state it was copied from and (2) the transition matrix of the demerged model
         */
        std::pair<std::vector<uint64_t>, storage::SparseMatrix<ValueType>> buildTransitionMatrix(automata::AcceptanceCondition const& acceptance, storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::SparseMatrix<ValueType> const& backwardTransitions);

        //std::unordered_map<std::string, RewardModelType> buildRewardModels(
        //storm::storage::SparseMatrix<ValueType> const& resultTransitionMatrix);
};
}
}
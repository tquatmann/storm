#pragma once

#include "DARewardProduct.h"
#include "storm/logic/Formula.h"
#include "storm/automata/acceptanceCondition.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/transformer/DAProductBuilder.h"

namespace storm {
namespace transformer {

template<typename ValueType, typename RewardModelType>
class DARewardProductBuilder {
    public:
        using Mdp = storm::models::sparse::Mdp<ValueType, RewardModelType>;

        DARewardProductBuilder(Mdp& productModel, std::vector<storm::automata::AcceptanceCondition> const& acceptanceConditions, Mdp const& originalModel): productModel(productModel), acceptanceConditions(acceptanceConditions), originalModel(originalModel) {}

        std::shared_ptr<DARewardProduct<ValueType>> build();

    private:
        using Row = typename storage::SparseMatrix<ValueType>::const_rows;
        Mdp& productModel;
        std::vector<storm::automata::AcceptanceCondition> const& acceptanceConditions;
        Mdp const& originalModel;
        uint64_t InvalidIndex = std::numeric_limits<uint64_t>::max();

        /*!
         * Adds a modified row to the matrix builder, where entries of states in the same MEC summed up
         * @param builder the matrix builder where the state-action pair is saved
         * @param stateToMec mapping of states to the MEC they are in or the numerical limit otherwise
         * @param row the state-action pair from the original model
         * @param numMecs the number of MECs in the original model
         */
        void modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, std::vector<uint64_t> const& stateToMec, Row row, uint64_t numMecs, bool builderEmpty=false);

        /*!
         * Processes the product by determining the quotient-MDP
         * @param transitionMatrix the transition matrix of the product
         * @param transitionMatrixBuilder the matrix builder for the modified matrix
         * @param stateToMec
         * @param mecs
         * @param choiceToModelChoice a mapping of choices in the modified model to ones in the original model
         * @return a mapping of each MEC to the actions that leave it with some probability
         */
        std::vector<std::list<uint64_t>> processProductMatrix(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, std::vector<uint64_t> const& stateToMec, storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs, std::vector<uint64_t>& choiceToModelChoice);

        /*!
         *
         * @param transitionMatrix
         * @param transitionMatrixBuilder
         * @param stateToMec
         * @param mecs
         * @param accEcs
         * @param mecsToLeavingActions
         * @return a list of actions that lead to an accepting end component in the modified model a.s.
         */
        std::list<uint64_t> addRepresentativeStates(storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                    storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, std::vector<uint64_t> const& stateToMec,
                                                    storage::MaximalEndComponentDecomposition<ValueType> const& mecs,
                                                    std::list<storage::MaximalEndComponent>& accEcs,
                                                    std::vector<std::list<uint64_t>> const& mecsToLeavingActions);

        /*!
         * Adds copies of the MACs to the modified model
         * @param transitionMatrix
         * @param transitionMatrixBuilder
         * @param stateToMec
         * @param mecs
         * @param accEcs
         */
        void addMACStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, uint64_t numberMECs, std::list<storage::MaximalEndComponent> accEcs);

        /*!
         *
         * @param transitionMatrix
         * @param mecs
         * @param accEcs
         */
        void computeConversionsFromModel(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs, std::list<storage::MaximalEndComponent> const& accEcs, std::vector<uint64_t>& stateToModelState, std::vector<uint64_t>& choiceToModelChoice);

        /*!
         * Computes the set of maximal accepting end components
         * @param acceptance the acceptance condition
         * @param transitionMatrix the transition matrix of the product
         * @param backwardTransitions the reversed transition relation
         * @return the list of maximal accepting end components
         */
        std::list<storage::MaximalEndComponent> computeAcceptingECs(automata::AcceptanceCondition const& acceptance,
                                                                    storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                    storm::storage::SparseMatrix<ValueType> const& backwardTransitions);


        static void removeSubsetEndComponents(std::list<storage::MaximalEndComponent>& endComponents);

        /*!
         * Lifts the initial states from the original model to the modified one
         * @param stateToMec mapping of states to the MEC they are in or the numerical limit otherwise
         * @param numberOfStates the number of states of the modified model
         * @return the initial states of the modified transition matrix
         */
        storm::storage::BitVector liftInitialStates(std::vector<uint64_t> const& stateToMec, uint64_t numberOfStates);
};
}
}
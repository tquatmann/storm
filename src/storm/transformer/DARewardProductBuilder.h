#pragma once

#include "DARewardProduct.h"
#include "storm/logic/Formula.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/transformer/DAProductBuilder.h"

namespace storm {
namespace transformer {

template<typename ValueType, typename RewardModelType>
class DARewardProductBuilder {
    public:
        using Mdp = storm::models::sparse::Mdp<ValueType, RewardModelType>;

        DARewardProductBuilder(DAProduct<Mdp>& product): product(product) {}

        std::shared_ptr<DARewardProduct<ValueType>> build();

    private:
        using Row = typename storage::SparseMatrix<ValueType>::rows;
        DAProduct<Mdp>& product;
        uint64_t InvalidIndex = std::numeric_limits<uint64_t>::max();

        /*!
         * Builds the transition matrix of the demerged product MDP.
         * @param acceptance the acceptance condition
         * @param transitionMatrix the transition matrix of the product
         * @param backwardTransitions the reversed transition relation
         * @return A pair consisting of (1) a map from each state in the new model to the state it was copied from and (2) the transition matrix of the demerged model
         */
        DARewardProduct<ValueType> buildTransitionMatrix(automata::AcceptanceCondition const& acceptance, storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::SparseMatrix<ValueType> const& backwardTransitions);

        /*!
         * Adds a modified row to the matrix builder, where entries of states in the same MEC summed up
         * @param builder the matrix builder where the state-action pair is saved
         * @param stateToMec mapping of states to the MEC they are in or the numerical limit otherwise
         * @param row the state-action pair from the original model
         * @param numMecs the number of MECs in the original model
         */
        void modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, std::vector<uint64_t>& stateToMec, Row row, uint64_t numMecs, bool builderEmpty=false);

        /*!
         * Computes the set of maximal accepting end components
         * @param acceptance the acceptance condition
         * @param transitionMatrix the transition matrix of the product
         * @param backwardTransitions the reversed transition relation
         * @return the list of maximal accepting end components
         */
        std::list<storage::MaximalEndComponent> computeAcceptingECs(automata::AcceptanceCondition const& acceptance,storm::storage::SparseMatrix<ValueType> const& transitionMatrix,storm::storage::SparseMatrix<ValueType> const& backwardTransitions);

};
}
}
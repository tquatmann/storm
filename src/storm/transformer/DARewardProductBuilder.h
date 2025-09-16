#pragma once

#include "DARewardProduct.h"
#include "storm/automata/AcceptanceCondition.h"
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

    DARewardProductBuilder(Mdp& productModel, std::vector<storm::automata::AcceptanceCondition::ptr> const& acceptanceConditions,
                           std::vector<uint64_t> const& productToModelState, Mdp const& originalModel)
        : productModel(productModel), acceptanceConditions(acceptanceConditions), productToModelState(productToModelState), originalModel(originalModel) {}

    std::shared_ptr<DARewardProduct<ValueType>> build();

   private:
    using Row = typename storage::SparseMatrix<ValueType>::const_rows;
    Mdp& productModel;
    std::vector<storm::automata::AcceptanceCondition::ptr> const& acceptanceConditions;
    Mdp const& originalModel;
    std::vector<uint64_t> const& productToModelState;
    uint64_t InvalidIndex = std::numeric_limits<uint64_t>::max();

    class MecDecompositionInfo {
       public:
        storm::storage::MaximalEndComponentDecomposition<ValueType> mecs;
        uint64_t numStatesInMecs;
        uint64_t numChoicesInMecs;
        std::vector<uint64_t> stateToMec;

        MecDecompositionInfo(storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs_, uint64_t numStatesInMecs_, uint64_t numChoicesInMecs_,
                             std::vector<uint64_t> const& stateToMec_)
            : mecs(mecs_), numStatesInMecs(numStatesInMecs_), numChoicesInMecs(numChoicesInMecs_), stateToMec(stateToMec_) {}
    };

    struct Conversions {
        std::vector<uint64_t> stateToModelState;
        std::vector<uint64_t> choiceToModelChoice;
        std::vector<uint64_t> modelStateToState;
    };

    /*!
     * Adds a modified row to the matrix builder, where entries of states in the same MEC summed up
     * @param builder the matrix builder where the state-action pair is saved
     * @param row the state-action pair from the original model
     * @param builderEmpty
     */
    void modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, Row row, MecDecompositionInfo const& mecDecompositionInfo,
                               Conversions& conversions, bool builderEmpty = false);

    /*!
     * Processes the product by determining the quotient-MDP
     * @param transitionMatrix the transition matrix of the product
     * @param transitionMatrixBuilder the matrix builder for the modified matrix
     * @return a mapping of each MEC to the actions that leave it with some probability
     */
    std::vector<std::list<uint64_t>> processProductMatrix(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                          storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
                                                          MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions);

    /*!
     *
     * @param transitionMatrix
     * @param transitionMatrixBuilder
     * @param accEcs
     * @param mecsToLeavingActions
     * @return a list of actions that lead to an accepting end component in the modified model a.s.
     */
    std::vector<std::list<uint64_t>> addRepresentativeStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                             storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
                                                             std::vector<std::list<storage::MaximalEndComponent>> const& accEcs,
                                                             std::vector<std::list<uint64_t>> const& mecsToLeavingActions,
                                                             MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions);

    /*!
     * Adds copies of the MACs to the modified model
     * @param transitionMatrix
     * @param transitionMatrixBuilder
     * @param accEcs
     */
    void addMACStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
                      std::list<storage::MaximalEndComponent> const& accEcs, MecDecompositionInfo mecDecompositionInfo);

    /*!
     *
     * @param transitionMatrix
     * @param accEcs
     */
    void computeConversionsFromModel(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                     std::vector<std::list<storage::MaximalEndComponent>> const& accEcs, MecDecompositionInfo const& mecDecompositionInfo,
                                     Conversions& conversions, uint64_t stateCounter, uint64_t choiceCounter);

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
    storm::storage::BitVector liftInitialStates(MecDecompositionInfo const& mecDecompositionInfo, uint64_t numberOfStates,
                                                storm::storage::BitVector const& initialStates) const;
};
}  // namespace transformer
}  // namespace storm
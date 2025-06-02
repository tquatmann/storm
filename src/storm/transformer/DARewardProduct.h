#pragma once

#include "storm/storage/SparseMatrix.h"

namespace storm {
namespace transformer {

template <typename ValueType>
class DARewardProduct {
    public:
        DARewardProduct(storage::SparseMatrix<ValueType> transitionMatrix, std::vector<uint64_t> const& stateToModelState, std::vector<uint64_t> const& actionToModelAction, std::vector<std::list<uint64_t>> const& reachingAccEcChoices, storm::storage::BitVector const& initialStates):
            transitionMatrix(transitionMatrix),
            stateToModelState(stateToModelState),
            actionToModelAction(actionToModelAction),
            reachingAccEcChoices(reachingAccEcChoices),
            initialStates(initialStates) {}

        storage::SparseMatrix<ValueType> getTransitionMatrix() {
            return transitionMatrix;
        }

        std::vector<uint64_t> getStateToModelState() {
            return stateToModelState;
        }

        storm::storage::BitVector getInitialStates() {
            return initialStates;
        }

        std::vector<uint64_t> getActionToModelAction() {
            return actionToModelAction;
        }

        std::vector<std::list<uint64_t>> getReachingAccEcChoices() {
            return reachingAccEcChoices;
        }

    private:
        storage::SparseMatrix<ValueType> transitionMatrix;
        std::vector<uint64_t> stateToModelState;
        std::vector<uint64_t> actionToModelAction;
        std::vector<std::list<uint64_t>> reachingAccEcChoices;
        storm::storage::BitVector initialStates;
};

}
}



#pragma once

#include "storm/storage/SparseMatrix.h"

namespace storm {
namespace transformer {

template <typename ValueType>
class DARewardProduct {
    public:
        DARewardProduct(storage::SparseMatrix<ValueType> transitionMatrix, std::vector<uint64_t> stateToModelState, std::vector<uint64_t> actionToModelAction, std::list<uint64_t> reachingAccEcChoices, storm::storage::BitVector initialStates):
            transitionMatrix(transitionMatrix),
            stateToModelState(stateToModelState),
            actionToModelAction(actionToModelAction),
            reachingAccEcChoices(reachingAccEcChoices) {}

        storage::SparseMatrix<ValueType> getTransitionMatrix() {
            return transitionMatrix;
        }

        std::vector<uint64_t> getStateToModelState() {
            return stateToModelState;
        }

    private:
        storage::SparseMatrix<ValueType> transitionMatrix;
        std::vector<uint64_t> stateToModelState;
        std::vector<uint64_t> actionToModelAction;
        std::list<uint64_t> reachingAccEcChoices;
        storm::storage::BitVector initialStates;
};

}
}



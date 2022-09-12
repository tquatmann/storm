#include "storm/environment/solver/MinMaxLpSolverEnvironment.h"

namespace storm {
void MinMaxLpSolverEnvironment::setUseEqualityForSingleActions(bool newValue) {
    useEqualityForSingleActions = newValue;
}
void MinMaxLpSolverEnvironment::setOptimizeOnlyForInitialState(bool newValue) {
    optimizeOnlyForInitialState = newValue;
}
void MinMaxLpSolverEnvironment::setUseNonTrivialBounds(bool newValue) {
    useNonTrivialBounds = newValue;
}

bool MinMaxLpSolverEnvironment::getUseEqualityForSingleActions() const {
    return useEqualityForSingleActions;
}
bool MinMaxLpSolverEnvironment::getOptimizeOnlyForInitialState() const {
    return optimizeOnlyForInitialState;
}
bool MinMaxLpSolverEnvironment::getUseNonTrivialBounds() const {
    return useNonTrivialBounds;
}
}  // namespace storm
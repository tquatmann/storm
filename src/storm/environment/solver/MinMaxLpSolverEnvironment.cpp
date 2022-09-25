#include "storm/environment/solver/MinMaxLpSolverEnvironment.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/MinMaxEquationSolverSettings.h"

namespace storm {

MinMaxLpSolverEnvironment::MinMaxLpSolverEnvironment() {
    auto const& minMaxSettings = storm::settings::getModule<storm::settings::modules::MinMaxEquationSolverSettings>();
    useNonTrivialBounds = minMaxSettings.getLpUseNonTrivialBounds();
    optimizeOnlyForInitialState = minMaxSettings.getLpUseOnlyInitialStateAsObjective();
    useEqualityForSingleActions = minMaxSettings.getLpUseEqualityForTrivialActions();
    optimizeFeasibilityQueries = minMaxSettings.getLpOptimizeFeasibilityQueries();
}

void MinMaxLpSolverEnvironment::setUseEqualityForSingleActions(bool newValue) {
    useEqualityForSingleActions = newValue;
}
void MinMaxLpSolverEnvironment::setOptimizeOnlyForInitialState(bool newValue) {
    optimizeOnlyForInitialState = newValue;
}
void MinMaxLpSolverEnvironment::setUseNonTrivialBounds(bool newValue) {
    useNonTrivialBounds = newValue;
}
void MinMaxLpSolverEnvironment::setOptimizeFeasibilityQueries(bool newValue) {
    optimizeFeasibilityQueries = newValue;
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
bool MinMaxLpSolverEnvironment::getOptimizeFeasibilityQueries() const {
    return optimizeFeasibilityQueries;
}
}  // namespace storm
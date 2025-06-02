#pragma once

#include <memory>

#include "storm/environment/Environment.h"
#include "storm/storage/BitVector.h"
#include "storm/logic/Formulas.h"
#include "storm/modelchecker/multiobjective/MultiObjectiveModelCheckingMethod.h"
#include "storm/modelchecker/results/CheckResult.h"

namespace storm {

class Environment;

namespace modelchecker {
namespace multiobjective {
typedef std::function<storm::storage::BitVector(storm::logic::Formula const&)> CheckFormulaCallback;

template<typename SparseModelType>
std::unique_ptr<CheckResult> performMultiObjectiveModelChecking(Environment const& env, SparseModelType const& model,
                                                                storm::logic::MultiObjectiveFormula const& formula);

}
}  // namespace modelchecker
}  // namespace storm

#pragma once

#include <memory>

#include "storm/environment/Environment.h"
#include "storm/logic/Formulas.h"
#include "storm/modelchecker/multiobjective/MultiObjectiveModelCheckingMethod.h"
#include "storm/modelchecker/results/CheckResult.h"
#include "storm/storage/BitVector.h"

namespace storm {

class Environment;

namespace modelchecker {
namespace multiobjective {
typedef std::function<storm::storage::BitVector(storm::logic::Formula const&)> CheckFormulaCallback;

template<typename SparseModelType>
std::unique_ptr<CheckResult> performMultiObjectiveModelChecking(Environment const& env, SparseModelType const& model,
                                                                storm::logic::MultiObjectiveFormula const& formula, bool produceScheduler = false);

}  // namespace multiobjective
}  // namespace modelchecker
}  // namespace storm

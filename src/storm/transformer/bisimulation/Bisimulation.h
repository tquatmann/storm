#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "storm/logic/FormulasForwardDeclarations.h"
#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Options.h"

namespace storm::bisimulation {

template<typename ValueType>
struct ReturnType {
    std::shared_ptr<storm::models::sparse::Model<ValueType>> quotient;
    std::vector<uint64_t> toQuotientStateMapping;
};

/*!
 * Computes the bisimulation quotient of the given model.
 *
 * @param model the model to minimize.
 * @param options the options to use for the minimization.
 * @param formulas the formulas that need to be preserved by the minimization.
 */
template<typename ValueType>
ReturnType<ValueType> applyBisimulationMinimization(storm::models::sparse::Model<ValueType> const& model, Options const& options = {},
                                                      std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas = {});

}  // namespace storm::bisimulation

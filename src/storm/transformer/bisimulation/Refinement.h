#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Signature.h"

namespace storm::bisimulation {

/*!
 * Performs splitter-based partition refinement.
 * @note only applicable to deterministic models.
 */
template<typename ValueType>
void performSplitterBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition, ValueType const tolerance);

/*!
 * Performs signature-based partition refinement.
 * @note only applicable to nondeterministic models.
 */
template<typename ValueType>
void performSignatureBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                     Signatures<ValueType>& signatures);

}  // namespace storm::bisimulation

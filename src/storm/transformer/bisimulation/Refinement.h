#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Signatures.h"

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
template<typename ValueType, SignatureMode Mode>
void performSignatureBasedRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                     Signatures<ValueType, Mode>& signatures);

}  // namespace storm::bisimulation

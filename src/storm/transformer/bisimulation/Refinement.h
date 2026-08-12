#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"

namespace storm::bisimulation {

template<typename ValueType>
void performPartitionRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition,
                                std::optional<std::vector<uint64_t>> const& choiceClasses = {});

}  // namespace storm::bisimulation

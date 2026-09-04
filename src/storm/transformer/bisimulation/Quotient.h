#pragma once

#include <memory>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Options.h"
#include "storm/transformer/bisimulation/PreservationInformation.h"
#include "storm/transformer/bisimulation/QuotientData.h"

namespace storm::bisimulation {

template<typename ValueType>
class Quotient {
   public:
    /*!
     * Builds the quotient model from the given partition (represented by quotientData), preserving the given labels/rewards.
     */
    static std::shared_ptr<storm::models::sparse::Model<ValueType>> buildFromPartition(
        storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Options const& options,
        storm::bisimulation::PreservationInformation const& preservationInformation, QuotientData<ValueType> const& quotientData);
};

}  // namespace storm::bisimulation

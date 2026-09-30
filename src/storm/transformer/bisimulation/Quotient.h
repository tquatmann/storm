#pragma once

#include <memory>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Options.h"
#include "storm/transformer/bisimulation/PreservationInformation.h"
#include "storm/transformer/bisimulation/QuotientData.h"

namespace storm::bisimulation {

template<typename ValueType, typename QuotientValueType = ValueType>
class Quotient {
   public:
    /*!
     * Builds the quotient model from the given partition (represented by quotientData), preserving the given labels/rewards.
     * @tparam ValueType the value type of the given model
     * @tparam QuotientValueType the type of the quotient model. Either coincides with ValueType or with IntervalType<ValueType> (for interval abstraction).
     * @note for weak bisimulation, quotientData must carry the weak bisimulation data and must have been built with the non-silent states as its preferred
     * representatives.
     */
    static std::shared_ptr<storm::models::sparse::Model<QuotientValueType>> buildFromPartition(
        storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Options const& options,
        storm::bisimulation::PreservationInformation const& preservationInformation, QuotientData<QuotientValueType> const& quotientData);
};

}  // namespace storm::bisimulation

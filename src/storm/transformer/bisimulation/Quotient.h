#pragma once

#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include "storm/adapters/IntervalForward.h"
#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Options.h"
#include "storm/transformer/bisimulation/PreservationInformation.h"

namespace storm::bisimulation {

class Partition;

enum class QuotientType { Exact, IntervalAbstraction };

template<QuotientType quotientType, typename ValueType>
class Quotient {
   public:
    using QuotientValueType = std::conditional_t<quotientType == QuotientType::IntervalAbstraction, IntervalType<ValueType>, ValueType>;

    /*!
     * The index mappings between an input model and its bisimulation quotient, derived from a partition of the input model's states.
     */
    struct IndexMapping {
        std::vector<uint64_t> toQuotientState;        // assigns to each input model state the corresponding quotient state
        std::vector<uint64_t> toRepresentativeState;  // assigns to each quotient state the corresponding representative state in the input model
        std::optional<std::vector<uint64_t>> toRepresentativeChoice;  // assigns to each quotient choice the corresponding representative choice in the input
                                                                      // model. Only present for nondeterministic models.
        // TODO needed if we deduplicate choices? std::optional<std::vector<uint64_t>> quotientStateToChoiceIndices;  // CSR-style mapping that specifies for
        // each quotient state the range of choice indices. Only present for nondeterministic models.
    };

    /*!
     * Computes the index mappings between the given model and the quotient described by the given partition.
     */
    static IndexMapping computeIndexMappings(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition const& partition,
                                             std::optional<std::vector<uint64_t>> const& choiceClasses = {});

    /*!
     * Builds the quotient model from the given partition (represented by indexMapping), preserving the given labels/rewards.
     */
    static std::shared_ptr<storm::models::sparse::Model<QuotientValueType>> buildFromPartition(
        storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Options const& options,
        storm::bisimulation::PreservationInformation const& preservationInformation, IndexMapping const& indexMapping);
};

}  // namespace storm::bisimulation

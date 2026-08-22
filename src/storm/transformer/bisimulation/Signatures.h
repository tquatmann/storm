#pragma once

#include <compare>
#include <cstdint>
#include <optional>
#include <set>
#include <vector>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Partition.h"

namespace storm::bisimulation {

enum class SignatureMode { Exact, Approximative };

/*!
 * The (static) signature of a state w.r.t. a partition, based on the distributions of its choices over blocks.
 * Supports both, exact and approximate signatures (for the latter, probabilities are equal up to a given tolerance).
 *  A state signature is seen as a set of choice signature.
 *  The set is ordered so that c_1 < c_2 < c_3 and c_1 ≈ c_3 implies c_1 ≈ c_2 ≈ c_3.
 */
template<typename ValueType_, SignatureMode Mode>
struct StateSignature {
    using ValueType = ValueType_;

    struct ChoiceSignature {
        uint64_t const choiceClass;                            // The class of this choice according to the provided choice classes.
        typename Partition::OrderedBlockMap<ValueType> distr;  // The accumulated transition values for each block.
        std::strong_ordering operator<=>(ChoiceSignature const& other) const;
        bool operator==(ChoiceSignature const& other) const;
        bool approxEqual(ChoiceSignature const& other, ValueType const tolerance) const
            requires(Mode == SignatureMode::Approximative);
    };
    std::map<ChoiceSignature, uint64_t> choices;  // Maps each choice to the choiceIndex it originated from
};

/*!
 * Computes and caches the state signatures
 */
template<typename ValueType, SignatureMode Mode>
class Signatures {
   public:
    using SignatureType = StateSignature<ValueType, Mode>;

    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition)
        requires(Mode == SignatureMode::Exact);

    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition, ValueType const& tolerance)
        requires(Mode == SignatureMode::Approximative);

    /*!
     * Updates the state signature of the given state index.
     * @param stateIndex
     */
    void updateStateSignature(uint64_t const stateIndex);

    struct SplitOrder {
        std::vector<SignatureType> const& signatures;
        bool operator()(uint64_t const state1, uint64_t const state2) const;
    };
    /*!
     * Returns a strict weak order over state indices based on their (most recently updated) signature.
     * @note the order might be out-of-date for states whose signature has not been updated since the most recent call to `updateStateSignature`
     */
    SplitOrder getSplitOrder() const;

    struct SplitCondition {
        std::vector<SignatureType> const& signatures;
        ValueType const tolerance;  // only meaningful for approximate signatures
        bool operator()(uint64_t const state1, uint64_t const state2) const;
    };

    /*!
     * @return a SplitCondition object sc that can be called on two states such that sc(state1, state2) returns true iff the two states should not be in the
     * same partition.
     */
    SplitCondition getSplitCondition() const;

   private:
    storm::models::sparse::Model<ValueType> const& model;
    Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;
    ValueType const tolerance;  // only meaningful for approximate signatures

    std::vector<SignatureType> stateSignatureCache;
};

}  // namespace storm::bisimulation

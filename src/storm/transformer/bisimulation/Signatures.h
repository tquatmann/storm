#pragma once

#include <compare>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/QuotientData.h"

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
        std::strong_ordering compare(ChoiceSignature const& other, [[maybe_unused]] ValueType const tolerance) const;
        bool operator==(ChoiceSignature const& other) const
            requires(Mode == SignatureMode::Exact);
        bool approxEqual(ChoiceSignature const& other, ValueType const tolerance) const
            requires(Mode == SignatureMode::Approximative);
    };
    std::vector<ChoiceSignature> choices;  // Maps each choice to the choiceIndex it originated from
};

/*!
 * Computes and caches the state signatures
 */
template<typename ValueType, SignatureMode Mode>
class Signatures {
   public:
    using StateSignature = StateSignature<ValueType, Mode>;
    using ChoiceSignature = typename StateSignature::ChoiceSignature;

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
        std::vector<StateSignature> const& signatures;
        ValueType const tolerance;  // only meaningful for approximate signatures
        bool operator()(uint64_t const state1, uint64_t const state2) const;
    };
    /*!
     * Returns a strict weak order over state indices based on their (most recently updated) signature.
     * @note the order might be out-of-date for states whose signature has not been updated since the most recent call to `updateStateSignature`
     */
    SplitOrder getSplitOrder() const;

    struct SplitCondition {
        std::vector<StateSignature> const& signatures;
        ValueType const tolerance;  // only meaningful for approximate signatures
        bool operator()(uint64_t const state1, uint64_t const state2) const;
    };

    /*!
     * @return a SplitCondition object sc that can be called on two states such that sc(state1, state2) returns true iff the two states should not be in the
     * same partition.
     */
    SplitCondition getSplitCondition() const;

    /*!
     * Fills in the choice mappings of the given (already state-mapped) quotient data: for every quotient state, the choices of its representative state
     * are deduplicated (exactly or up to tolerance, depending on Mode) into the quotient choices, and every original model choice is mapped to the
     * quotient choice that represents it.
     * @note only applicable to nondeterministic models.
     */
    void extendQuotientData(QuotientData<ValueType>& quotientData) const;

   private:
    struct ChoiceOrder {
        ValueType const tolerance;  // only meaningful for approximate signatures
        bool operator()(ChoiceSignature const& choice1, ChoiceSignature const& choice2) const;
    };

    /*!
     * Merges the choices of the given state by their signature: choices whose signature compares equal (Exact mode) or approximately equal within
     * tolerance (Approximative mode) are merged into the same group. Used by both updateStateSignature (which only keeps one representative choice index
     * per group) and moveToQuotientData (which needs every original choice index per group).
     */
    std::map<ChoiceSignature, std::vector<uint64_t>, ChoiceOrder> mergeEquivalentChoices(uint64_t const stateIndex) const;

    storm::models::sparse::Model<ValueType> const& model;
    Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;
    ValueType const tolerance;  // only meaningful for approximate signatures

    std::vector<StateSignature> stateSignatureCache;
};

}  // namespace storm::bisimulation

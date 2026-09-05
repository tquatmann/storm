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
        uint64_t choiceClass;                                  // The class of this choice according to the provided choice classes.
        typename Partition::OrderedBlockMap<ValueType> distr;  // The accumulated transition values for each block.
        std::strong_ordering compare(ChoiceSignature const& other, [[maybe_unused]] ValueType const tolerance) const;
        bool operator==(ChoiceSignature const& other) const
            requires(Mode == SignatureMode::Exact);
        bool approxEqual(ChoiceSignature const& other, ValueType const tolerance) const
            requires(Mode == SignatureMode::Approximative);
    };

    /*!
     * Find the given choice signature in the list.
     * Returns an (iterator,bool)-pair. The bool is true iff the given choice signature is contained.
     * The iterator is either (exact/approximately) equal to the given signature, or the position where it would be inserted.
     * @return
     */
    auto find(ChoiceSignature const& signature, ValueType const tolerance) const {
        auto it =
            std::lower_bound(choices.begin(), choices.end(), signature, [this, &tolerance](ChoiceSignature const& choice1, ChoiceSignature const& choice2) {
                return choice1.compare(choice2, tolerance) == std::strong_ordering::less;
            });
        if constexpr (Mode == SignatureMode::Exact) {
            return std::make_pair(it, it != choices.end() && *it == signature);
        } else {
            // For approximate signatures, we need to check if the signature is approximately equal to any neighbouring signature.
            if (it != choices.end() && it->approxEqual(signature, tolerance)) {
                return std::make_pair(it, true);
            }
            if (it != choices.begin() && std::prev(it)->approxEqual(signature, tolerance)) {
                return std::make_pair(std::prev(it), true);
            }
            return std::make_pair(it, false);
        }
    }

    std::vector<ChoiceSignature> choices;  // Sorted
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
    /*!
     * Gets the signature of the given choice index.
     * The underlying distribution maps blocks to (exact) aggregated probability.
     */
    ChoiceSignature getChoiceSignature(uint64_t const choiceIndex) const;

    storm::models::sparse::Model<ValueType> const& model;
    Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;

    /*!
     * Half of the tolerance passed to the constructor; only meaningful for approximate signatures.
     *
     * This is used to distribute the tolerance from the constructor across two approximate-equality (≈) checks.
     * To see why taking half of the tolerance is necessary, consider the following example:
     * Let s_1 and s_2 be two states of a block. Let c_1 be a choice of s_1 and let c_2 and c'_2 be two choices of s_2 such that c_1 ≈ c_2 ≈ c'_2.
     * Then, c'_2 does not need to appear in s_2's signature as it is represented by c_2.
     * The quotient might only have s_1 and c_1 as representatives for the block but c_1 ≈ c'_2 does not need to hold. c_1 ≈ c_2 ≈ c'_2 only imples that
     * c_1 and c'_2 are within twice of the tolerance considered for ≈.
     */
    ValueType const halfTolerance;

    std::vector<StateSignature> stateSignatureCache;
};

}  // namespace storm::bisimulation

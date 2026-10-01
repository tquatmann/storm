#pragma once

#include <compare>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/QuotientData.h"

namespace storm::bisimulation {

enum class SignatureMode {
    Exact,                 // Two values are equal iff they coincide.
    Approximative,         // Two values are equal iff they differ by at most a given tolerance.
    IntervalAbstraction  // The signature values are intervals that cover all values of the represented states.
};

/*!
 * Computes and caches the state signatures and provides state-based order / condition for splitting in signature refinement.
 * @tparam ValueType the type of the transition values in the model
 * @tparam Mode the mode of the signature computation, cf. SignatureMode
 * @tparam QuotientValueType the type of the quotient model. Either coincides with ValueType or with IntervalType<ValueType> (for interval abstraction).
 */
template<typename ValueType, SignatureMode Mode, typename QuotientValueType = ValueType>
class Signatures {
    static_assert(Mode != SignatureMode::IntervalAbstraction || storm::IsIntervalType<QuotientValueType>,
                  "SignatureMode::IntervalAbstraction requires an interval quotient value type.");
    static_assert(Mode != SignatureMode::Approximative || !storm::IsIntervalType<QuotientValueType>,
                  "Approximative signatures over interval values are expressed by SignatureMode::IntervalAbstraction.");

   public:
    using ToleranceType = storm::IntervalBaseType<QuotientValueType>;
    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition)
        requires(Mode == SignatureMode::Exact);

    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition, ToleranceType const& tolerance)
        requires(Mode != SignatureMode::Exact);

    /*!
     * Updates the state signature of the given state index.
     * @param stateIndex
     */
    void updateStateSignature(uint64_t const stateIndex);

   private: // Forward-declared structs; defined privately further down.
    struct StateSignature;
    struct TemporaryStateSignature;

   public:
    /*!
     * Represents a strict weak order over state indices based on their signature.
     * @note assumes that only states are compared for which updateStateSignature has been called since the last change of the partition.
     */
    class SplitOrder {
       public:
        explicit SplitOrder(std::vector<StateSignature> const& signatures) : signatures(signatures) {}
        bool operator()(uint64_t const state1, uint64_t const state2) const;

       private:
        std::vector<StateSignature> const& signatures;
    };

    /*!
     * Instantiates a SplitOrder whose ties exactly coincide with true signature equality (operator==)
     */
    SplitOrder getEquivalenceSplitOrder() const
        requires(Mode == SignatureMode::Exact);

    /*!
     * Instantiates a SplitOrder that only compares with respect to structural equality of the state signature (ChoiceSignature::compareStructure)
     */
    SplitOrder getStructuralSplitOrder() const
        requires(Mode != SignatureMode::Exact);

    /*!
     * Represents a condition for splitting two states based on their signature.
     * @note assumes that only states are compared for which updateStateSignature has been called since the last change of the partition.
     * @note assumes that the two states are equivalent with respect to getStructuralSplitOrder()
     * operator() returns true iff the anchor state and the candidate state should *not* be grouped together.
     */
    struct SplitCondition {
        SplitCondition(std::vector<StateSignature>& signatures, ToleranceType const& tolerance, TemporaryStateSignature& tmpStateSignature)
            : signatures(signatures), tolerance(tolerance), tmpStateSignature(tmpStateSignature) {}
        bool operator()(uint64_t const anchorState, uint64_t const candidateState) const
            requires(Mode == SignatureMode::Approximative);
        bool operator()(uint64_t const anchorState, uint64_t const candidateState)
            requires(Mode == SignatureMode::IntervalAbstraction);

       private:
        std::vector<StateSignature>& signatures;
        ToleranceType const tolerance;
        TemporaryStateSignature& tmpStateSignature;
    };

    /*!
     * Instantiates a SplitCondition.
     */
    SplitCondition getClusteringSplitCondition()
        requires(Mode != SignatureMode::Exact);

    /*!
     * Copies the state signature from src to dst. Assumes that both have the same choice count and the choices are pairwise structurally equal (compareStructure() == 0).
     * @note: currently only needed for IntervalAbstraction.
     */
    void copyStructuralEquivalentStateSignature(uint64_t const& srcState, uint64_t const& dstState) requires(Mode == SignatureMode::IntervalAbstraction);

    /*!
     * Fills in the choice mappings of the given (already state-mapped) quotient data: for every quotient state, the choices of its representative state
     * are deduplicated (exactly or up to tolerance, depending on Mode) into the quotient choices.
     * @param createQuotientChoiceMapping if set, every original model choice is additionally mapped to the quotient choice that represents it
     *   (SignatureData::toQuotientChoice). For each original model state, this mapping is surjective onto the choices of the corresponding quotient state,
     *   so that, e.g., a scheduler for the quotient can be translated back to the original model.
     *   If not set, that (comparatively expensive) mapping remains empty.
     * @throws UnexpectedException if the choices of a state cannot be related to the choices of the state representing it in the quotient, which can happen
     *   in approximative mode since approximate equality is not transitive.
     */
    void extendQuotientData(QuotientData<QuotientValueType>& quotientData, bool const createQuotientChoiceMapping) const;

   private:
    // For all choices, we store the current distributions over successor blocks in choiceDistributionStorage and use a view into that for ChoiceSignatures.
    using BlockDistributionView = std::span<std::pair<Partition::Block, QuotientValueType>>;

    /*!
     * The signature of a single choice, consisting of its choice class and the distribution over successor blocks (Partition::Block -> QuotientValueType).
     */
    struct ConcreteChoiceSignature {
        uint64_t choiceClass{0};  // The class of this choice according to the provided choice classes.
        // The accumulated transition values for each block, sorted by Partition::BlockCompare. This is a (non-owning) view into the choice's permanent
        // slot of Signatures::choiceDistributionStorage, see there.
        BlockDistributionView distr;

        // Orders by choiceClass and the set of blocks (identity and order) of distr; independent of the actual transition values.
        std::strong_ordering compareStructure(ConcreteChoiceSignature const& other) const;
        // strong_ordering in exact mode (equal signatures are substitutable); only weak_ordering in approximative/interval abstraction mode, since ties don't
        // imply approxEqual.
        using ComparisonResult = std::conditional_t<Mode == SignatureMode::Exact, std::strong_ordering, std::weak_ordering>;
        // Strict weak order for sorting/lower_bound. In approximative/interval abstraction mode, ties don't imply (approximate) equality - see
        // StateSignature::find.
        ComparisonResult compare(ConcreteChoiceSignature const& other) const;
        // True iff the two signatures have the same structure and their distr values pairwise differ by at most the given tolerance.
        bool approximatelyEqual(ConcreteChoiceSignature const& other, ToleranceType const& tolerance) const
            requires(Mode == SignatureMode::Approximative);
    };
    struct AbstractChoiceSignature : public ConcreteChoiceSignature {
        // The lower and upper delta values accumulate the amount that the choice signature got widened compared to the concrete signature
        ToleranceType lowerDelta{storm::utility::zero<ToleranceType>()}, upperDelta{storm::utility::zero<ToleranceType>()};

        /*!
         * Widens the current choice signature to also include the other choice signature.
         * The lowerDelta and upperDelta values keep track of how much the intervals were widened compared to the concrete signatures.
         *
         * @param other the other choice signature. Assumes that the two signatures have the same structure (compareStructure() == 0).
         */
        void enhance(AbstractChoiceSignature const& other) requires(Mode == SignatureMode::IntervalAbstraction);

        /*!
         * We call two choice signatures c1 and c2 compatible if the delta values of c1.enhance(c2) are both less than or equal to the tolerance.
         * @param other the other choice signature. Assumes that the two signatures have the same structure (compareStructure() == 0).
         * @param requireContainsOther if true, the other signature must additionally be contained in the current signature in order to be compatible.
         */
        bool isCompatibleWith(AbstractChoiceSignature const& other, ToleranceType const& tolerance, bool requireContainsOther) const
            requires(Mode == SignatureMode::IntervalAbstraction);
    };
    using ChoiceSignature = std::conditional_t<Mode == SignatureMode::IntervalAbstraction, AbstractChoiceSignature, ConcreteChoiceSignature>;

    using ChoiceSignatureIterator = std::vector<ChoiceSignature>::const_iterator;

    /*!
     * The (static) signature of a state w.r.t. a partition.
     *  A state signature is seen as an ordered set of deduplicated choice signatures, where the order is defined by ChoiceSignature::compare and the
     *  deduplication is defined by the equality used for the mode (exact equality for Exact, approximate equality otherwise, see ChoiceSignature).
     */
    struct StateSignature {
        /*!
         * Find the given choice signature in the list.
         * Returns an (iterator,bool)-pair. The bool is true iff the given choice signature is contained.
         * The iterator is either (exact/approximately) equal to the given signature, or the position where it would be inserted.
         * Two choice signatures are approximately equal (w.r.t. a tolerance) iff they have equal structure (compareStructure) and their distr values pairwise
         * differ by at most the tolerance.
         * @param requireContainsSignature only relevant for Mode==IntervalAbstraction: if true, the found signature must additionally contained the given signature.
         */
        std::pair<ChoiceSignatureIterator, bool> find(ChoiceSignature const& signature, ToleranceType const& tolerance, [[maybe_unused]] bool requireContainsSignature = false)  const;

        /*!
         * Like find(), but starts searching at `hint` instead of locating the window of candidates via lower_bound. The window consists of the entries that are
         * compareStructure-equal to signature and whose first entry is within tolerance of the first entry of signature.
         * @param hint any iterator into choices, including choices.end(). The search moves linearly from hint to the window, so it is fast if hint is close to
         * it.
         * @note requires !signature.distr.empty().
         * @note returns {hint, false} if signature is not found.
         *  @param requireContainsSignature only relevant for Mode==IntervalAbstraction: if true, the found signature must additionally contained the given signature.
         */
        std::pair<ChoiceSignatureIterator, bool> findWithHint(ChoiceSignatureIterator const hint, ChoiceSignature const& signature,
                                                              ToleranceType const& tolerance,  [[maybe_unused]] bool requireContainsSignature = false) const
            requires(Mode != SignatureMode::Exact);

        void insert(ChoiceSignature const& choiceSignature, ToleranceType const& tolerance);

        std::vector<ChoiceSignature> choices;  // Always ordered and deduplicated as described above.
    };

    /*!
     * A small helper class that holds a state signature as a deep copy.
     * Currently only used for IntervalAbstraction signatures.
     */
    class TemporaryStateSignature {
    public:
        TemporaryStateSignature(uint64_t const maxTransitions) : distrStorage(maxTransitions) {}

        StateSignature& load(StateSignature const& src) requires(Mode == SignatureMode::IntervalAbstraction) {
            uint64_t offset = 0;
            sig.choices.resize(src.distr.size());
            for (uint64_t choiceIndex = 0; choiceIndex < src.choices.size(); ++choiceIndex) {
                auto& sigChoice = sig.choices[choiceIndex];
                auto const& srcChoice = src.choices[choiceIndex];
                sigChoice.choiceClass = srcChoice.choiceClass;
                sigChoice.lowerDelta = srcChoice.lowerDelta;
                sigChoice.upperDelta = srcChoice.upperDelta;
                sigChoice.distr = BlockDistributionView(distrStorage.data() + offset, srcChoice.distr.size()); // Initialize the distr view
                std::copy(srcChoice.distr.begin(), srcChoice.distr.end(), sigChoice.distr.begin()); // Copy the distribution
                offset += srcChoice.distr.size();
            }
            return sig;
        }

        void store(StateSignature& dst) {
            STORM_LOG_ASSERT(sig.choices.size() == dst.distr.size(), "Expected that the destination has the same choice count.");
            for (uint64_t choiceIndex = 0; choiceIndex < sig.choices.size; ++choiceIndex) {
                auto const& sigChoice = sig.choices[choiceIndex];
                auto& dstChoice = dst.choices[choiceIndex];
                dstChoice.choiceClass = sigChoice.choiceClass;
                dstChoice.lowerDelta = sigChoice.lowerDelta;
                dstChoice.upperDelta = sigChoice.upperDelta;
                STORM_LOG_ASSERT(sigChoice.distr.size() == dstChoice.distr.size(), "Expected that the destination has the same distribution size.");
                std::copy(sigChoice.distr.begin(), sigChoice.distr.end(), dstChoice.distr.begin());
            }
        }

    private:
        StateSignature sig;
        std::vector<std::pair<Partition::Block, QuotientValueType>> distrStorage;
    };

    /*!
     * Assigns to every choice of the given state the quotient choice that represents it, i.e., fills the corresponding entries of toQuotientChoice.
     * The assignment is surjective on the quotient choices of the corresponding quotient state, so that, e.g., a scheduler for the quotient can be translated
     * back to this state. Moreover, the values of a choice and those of the quotient choice representing it differ by at most the tolerance this instance was
     * created with.
     * @throws UnexpectedException if no such assignment is found, cf. extendQuotientData.
     * @param state a state of the block that representativeSignature belongs to
     * @param representativeSignature the signature of the representative state of that block
     * @param choiceSignatureToQuotientChoiceIndex assigns to each index of representativeSignature.choices the quotient choice derived from it
     * @param toQuotientChoice out-parameter, cf. QuotientData::toQuotientChoice
     */
    void addQuotientChoiceMapping(uint64_t const state, StateSignature const& representativeSignature,
                                  std::vector<uint64_t> const& choiceSignatureToQuotientChoiceIndex, std::vector<uint64_t>& toQuotientChoice) const;

    /*!
     * Builds the signature of the given choice index.
     * The underlying distribution maps blocks to aggregated probability.
     * @note the returned ChoiceSignature::distr is a view into choiceDistributionStorage; calling buildChoiceSignature(choiceIndex) again later overwrites it.
     */
    ChoiceSignature buildChoiceSignature(uint64_t const choiceIndex);

    /*!
     * Gets the signature of the given choice index. The signature must have been built before.
     * @note the returned ChoiceSignature::distr is a view into choiceDistributionStorage; calling buildChoiceSignature(choiceIndex) later overwrites it.
     */
    ChoiceSignature getChoiceSignature(uint64_t const choiceIndex) const;

    /*!
     * A reusable, sparse accumulator for building the distribution (Partition::Block -> ValueType) of a single choice, without allocating a fresh map
     * for every choice: values is a dense array indexed by a block's front() element (unique per block, since blocks are disjoint), and support records which
     * blocks were touched so that extract() only has to look at (and reset) those, not the whole array.
     * @note the cache must not be used across a change of the partition, since a block's front() element is only a stable identifier for that block as long
     * as the partition does not change.
     */
    class ChoiceSignatureCache {
       public:
        explicit ChoiceSignatureCache(uint64_t const numStates);

        /*!
         * Adds value to the accumulated value of the block b. Registers b as touched the first time this happens.
         */
        void addValue(Partition::Block const& block, ValueType const& value);

        /*!
         * Writes the accumulated (block, value) pairs, sorted by Partition::BlockCompare, into the given range (which must be large enough to hold them)
         * @return a sub-view over the entries actually written to dest.
         * @note resets the cache (i.e. all accumulated values) so that it can be reused for the next choice.
         */
        BlockDistributionView extract(BlockDistributionView const dest);

       private:
        std::vector<ValueType> values;
        std::vector<Partition::Block> support;
    } choiceSignatureCache;

    storm::models::sparse::Model<ValueType> const& model;
    Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;

    /*!
     * Backing storage for all ChoiceSignature::distr views. Choice i is permanently assigned the sub-range starting at the same offset as row i within
     * model's transition matrix (i.e. distance(model.getTransitionMatrix().begin(), model.getTransitionMatrix().begin(i))), which is guaranteed to fit
     * distr's entries since distr can have at most as many entries as row i of the transition matrix (deduplicating by block only ever shrinks it). This
     * avoids allocating a fresh container for every buildChoiceSignature call. Sized to fit every row's entries at once, since (row-)offsets are permanent.
     */
    mutable std::vector<std::pair<Partition::Block, QuotientValueType>> choiceDistributionStorage;

    /*!
     * Half of the tolerance passed to the constructor; only meaningful for approximate signatures.
     *
     * This is used to distribute the tolerance from the constructor across two approximate-equality (≈) checks.
     * To see why taking half of the tolerance is necessary, consider the following example:
     * Let s_1 and s_2 be two states of a block. Let c_1 be a choice of s_1 and let c_2 and c'_2 be two choices of s_2 such that c_1 ≈ c_2 ≈ c'_2.
     * Then, c'_2 does not need to appear in s_2's signature as it is represented by c_2.
     * The quotient might only have s_1 and c_1 as representatives for the block but c_1 ≈ c'_2 does not need to hold. c_1 ≈ c_2 ≈ c'_2 only implies that
     * c_1 and c'_2 are within twice of the tolerance considered for ≈.
     */
    ToleranceType const halfTolerance;

    std::vector<StateSignature> stateSignatureCache;

    /*!
     * Scratch space for a temporary state signature.
     */
    TemporaryStateSignature tmpStateSignature;
};

}  // namespace storm::bisimulation

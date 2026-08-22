#pragma once

#include <cstdint>
#include <optional>
#include <set>
#include <vector>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Partition.h"

namespace storm::bisimulation {

/*!
 * The (static) signature of a state w.r.t. a partition, based on the exact distributions of its choices.
 */
template<typename ValueType_>
struct ExactStateSignature {
    using ValueType = ValueType_;
    static constexpr bool IsExact = true;

    struct ChoiceSignature {
        uint64_t const choiceClass;
        typename Partition::OrderedBlockMap<ValueType> distr;
        bool operator<(ChoiceSignature const& other) const;
    };
    std::set<ChoiceSignature> choices;
};

template<typename StateSignature>
class Signatures;

/*!
 * A strict weak order over state indices, based on their (most recently updated) signature.
 */
template<typename StateSignature>
struct SplitOrder {
    Signatures<StateSignature> const& signatures;
    bool operator()(uint64_t const state1, uint64_t const state2) const;
};

/*!
 * Computes and caches the state signatures
 */
template<typename StateSignature>
class Signatures {
   public:
    using ValueType = typename StateSignature::ValueType;
    static constexpr bool IsExact = StateSignature::IsExact;

    friend struct SplitOrder<StateSignature>;

    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition);

    /*!
     * Updates the state signature of the given state index.
     * @param stateIndex
     */
    void updateStateSignature(uint64_t const stateIndex);

    /*!
     * Returns a strict weak order over state indices based on their (most recently updated) signature.
     * @note the order might be out-of-date for states whose signature has not been updated since the most recent call to `updateStateSignature`
     * @note only available if IsExact is true, since only exact signatures are guaranteed to induce a strict weak order.
     */
    SplitOrder<StateSignature> getExactSplitOrder() const
        requires IsExact;

   private:
    storm::models::sparse::Model<ValueType> const& model;
    Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;

    std::vector<StateSignature> stateSignatureCache;
};

}  // namespace storm::bisimulation

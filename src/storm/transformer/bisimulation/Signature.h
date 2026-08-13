#pragma once

#include <cstdint>
#include <optional>
#include <set>
#include <vector>

#include "storm/models/sparse/ModelForward.h"
#include "storm/transformer/bisimulation/Partition.h"

namespace storm::bisimulation {

/*!
 * Computes and caches the state signatures
 */
template<typename ValueType>
class Signatures {
   public:
/*!
     * The signature of a choice w.r.t. a partition: the class the choice belongs to together with the distribution over the blocks of the partition.
     */
    struct ChoiceSignature {
        uint64_t const choiceClass;
        Partition::OrderedBlockMap<ValueType> distr;  // todo: try alternatives

        // TODO: handle tolerance!
        bool operator<(ChoiceSignature const& other) const;
        bool operator==(ChoiceSignature const& other) const;
    };

    /*!
     * The (static) signature of a state w.r.t. a partition: the set of signatures of its choices.
     */
    struct StateSignature {
        std::set<ChoiceSignature> choices;

        bool operator<(StateSignature const& other) const;
    };

    Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
               storm::bisimulation::Partition const& partition);

    /*!
     * Updates the state signature of the given state index.
     * @param stateIndex
     */
    void updateStateSignature(uint64_t const stateIndex);

    /*!
     * Returns the signature of a state w.r.t. the partition.
     * @note the returned signature might be out-of-date if the partition has been modified since the most recent call to `updateStateSignature`
     */
    StateSignature const& getStateSignature(uint64_t const stateIndex) const;

   private:
    /*!
     * Computes the choice signature for the given choice index.
     * @param choiceIndex
     * @return
     */
    ChoiceSignature computeChoiceSignature(uint64_t const choiceIndex) const;

    storm::models::sparse::Model<ValueType> const& model;
    storm::bisimulation::Partition const& partition;
    std::optional<std::vector<uint64_t>> const& choiceClasses;

    std::vector<StateSignature> stateSignatureCache;
};

}  // namespace storm::bisimulation

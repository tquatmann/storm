#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <vector>

#include "storm/logic/FormulasForwardDeclarations.h"
#include "storm/models/sparse/ModelForward.h"
#include "storm/storage/BitVector.h"
#include "storm/transformer/bisimulation/Options.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/PreservationInformation.h"

namespace storm::bisimulation {

template<typename ValueType>
class Initialization {
   public:
    Initialization(storm::models::sparse::Model<ValueType> const& model, Options const& options,
                   std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);

    PreservationInformation getPreservationInformation() const;

    std::optional<std::vector<uint64_t>> getChoiceClasses() const;

    Partition getInitialStatePartition(std::optional<std::vector<uint64_t>> const& choiceClasses = {}) const;

   private:
    storm::models::sparse::Model<ValueType> const& model;
    Options const options;
    std::vector<std::shared_ptr<storm::logic::Formula const>> const formulas;

    /*!
     * Bookkeeping of the labels/rewards that need to be preserved by the bisimulation, split into their boolean-, integer-, and value-typed parts.
     */
    struct PreservedAnnotations {
        std::vector<std::reference_wrapper<storm::storage::BitVector const>> booleans;
        std::vector<std::span<uint64_t const>> integers;
        std::vector<std::span<ValueType const>> values;

        /*!
         * @return true iff there are no preserved annotations
         */
        bool empty() const;

        /*!
         * Splits all blocks in the partition with respect to the stored annotations.
         * @post Each two elements in a block of the partition have the same annotations.
         * @param partition The partition to be refined
         * @param extraAnnotation if non-empty, the partition is also split according to this additional annotation
         */
        void applySplit(Partition& partition, std::vector<uint64_t> const& extraAnnotation = {}) const;
    } preservedStateAnnotations, preservedChoiceAnnotations;
};

}  // namespace storm::bisimulation

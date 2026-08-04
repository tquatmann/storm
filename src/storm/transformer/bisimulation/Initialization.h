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
#include "storm/transformer/bisimulation/PreservationInformation.h"

namespace storm::bisimulation {

template<typename ValueType>
class Initialization {
   public:
    Initialization(storm::models::sparse::Model<ValueType> const& model, Options const& options,
                   std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);

    PreservationInformation getPreservationInformation() const;

    std::vector<uint64_t> getChoiceClasses() const;

    std::vector<uint64_t> getStateClasses(std::vector<uint64_t> const& choiceClasses) const;

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
    };
    PreservedAnnotations preservedStateAnnotations;
    PreservedAnnotations preservedChoiceAnnotations;

    /*!
     * Compares two elements (states or choices) by their preserved annotations, using extraAnnotation (e.g. a local choice index, or a state's ordered
     * choice classes) as an additional, higher-priority criterion.
     */
    struct LocalSignatureComp {
        LocalSignatureComp(PreservedAnnotations const& annotations, std::vector<uint64_t> const& extraAnnotation = {})
            : annotations(annotations), extraAnnotation(extraAnnotation) {}

        bool operator()(uint64_t a, uint64_t b) const;

        PreservedAnnotations const& annotations;
        std::vector<uint64_t> const& extraAnnotation;
    };
};

}  // namespace storm::bisimulation

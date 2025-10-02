#pragma once

#include <storm/storage/MaximalEndComponentDecomposition.h>

#include "storm/automata/Automaton.h"
#include "storm/storage/BitVector.h"
#include "storm/transformer/DAProduct.h"
#include "storm/transformer/DAProductBuilder.h"
#include "storm/transformer/Product.h"
#include "storm/transformer/ProductBuilder.h"

#include <vector>

namespace storm::logic {
class PathFormula;
}
namespace storm {
namespace transformer {
class LDBAProductBuilder {
   public:
    LDBAProductBuilder(const storm::automata::LimitDeterministicAutomaton& ldba, const std::vector<storm::storage::BitVector>& statesForAP)
        : ldba(ldba), statesForAP(statesForAP) {}

    template<typename ValueType>
    typename DAProduct<ValueType>::ptr build(storm::models::sparse::Model<ValueType> const& originalModel) const {
        // First build the product model which will be equipped with the initial states but no other label
        auto result = build(originalModel.getTransitionMatrix(), originalModel.getType(), originalModel.getInitialStates(), "init");
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                        "LDBA products currently do not support lifting original model information to the product.");
        // TODO: the LDBA product works a bit differently compared to the DA product because the current automaton state is based on a nondeterministic choice.
        // need to adapt this
        STORM_LOG_THROW(!originalModel.isOfType(storm::models::ModelType::MarkovAutomaton), storm::exceptions::NotSupportedException,
                        "LDBA products with Markov automata are currently not supported as they might introduce nondeterministic choices at Markovian states.");

        DAProductBuilder::liftOriginalModelInformationToProduct(originalModel, result);
        return result;
    }

    template<typename ValueType>
    typename DAProduct<ValueType>::ptr build(const storm::storage::SparseMatrix<ValueType>& transitionMatrix, storm::models::ModelType const& modelType,
                                             const storm::storage::BitVector& statesOfInterest, std::string const& statesOfInterestLabel = "soi") const {
        storage::MaximalEndComponentDecomposition mecs(transitionMatrix, transitionMatrix.transpose(true));
        storm::storage::BitVector mecStates(transitionMatrix.getRowGroupCount());
        for (auto const& ec : mecs) {
            for (auto const& [state, _] : ec) {
                mecStates.set(state, true);
            }
        }

        /*
        std::cout << "Automaton" << std::endl;
        for (int i = 0; i < ldba.getNumberOfStates(); i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << "(" << i << "," << j << "): " << ldba.getSuccessors(i, j) << std::endl;
            }
        }

        std::cout << "Labels" << std::endl;
        for (int i = 0; i < transitionMatrix.getRowGroupCount(); i++) {
            std::cout << i << ": " << getLabelForState(i) << std::endl;
        }
        */
        typename Product<ValueType>::ptr product =
            ProductBuilder<ValueType>::buildLimitDeterministicProduct(transitionMatrix, modelType, *this, statesOfInterest, statesOfInterestLabel, mecStates);
        storm::automata::AcceptanceCondition::ptr prodAcceptance = ldba.getAcceptance()->lift(
            product->getProductModel().getNumberOfStates(), [&product](std::size_t prodState) { return product->getAutomatonState(prodState); });

        return typename DAProduct<ValueType>::ptr(new DAProduct<ValueType>(std::move(*product), prodAcceptance));
    }

    storm::storage::sparse::state_type getInitialState(storm::storage::sparse::state_type MdpState) const {
        return ldba.getInitialState();
    }

    const std::vector<storm::storage::sparse::state_type> getSuccessors(storm::storage::sparse::state_type automatonFrom,
                                                                        storm::storage::sparse::state_type modelTo) const {
        std::vector<storm::storage::sparse::state_type> succs;
        for (auto const& succ : ldba.getSuccessors(automatonFrom, getLabelForState(modelTo))) {
            succs.push_back(succ);
        }
        return succs;
    }

    const storm::storage::BitVector getAcceptingPart() const {
        return ldba.computeAcceptingPart();
    }

   private:
    const storm::automata::LimitDeterministicAutomaton& ldba;
    const std::vector<storm::storage::BitVector>& statesForAP;

    storm::automata::APSet::alphabet_element getLabelForState(storm::storage::sparse::state_type s) const {
        storm::automata::APSet::alphabet_element label = ldba.getAPSet().elementAllFalse();
        for (unsigned int ap = 0; ap < ldba.getAPSet().size(); ap++) {
            if (statesForAP.at(ap).get(s)) {
                label = ldba.getAPSet().elementAddAP(label, ap);
            }
        }
        return label;
    }
};
}  // namespace transformer
}  // namespace storm

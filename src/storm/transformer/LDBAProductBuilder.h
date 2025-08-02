#pragma once

#include <storm/storage/MaximalEndComponentDecomposition.h>

#include "storm/automata/DeterministicAutomaton.h"
#include "storm/storage/BitVector.h"
#include "storm/transformer/DAProduct.h"
#include "storm/transformer/Product.h"
#include "storm/transformer/ProductBuilder.h"

#include <vector>

#include "storm/automata/limitdeterministicautomata/LimitDeterministicAutomaton.h"

namespace storm::logic {
class PathFormula;
}
namespace storm {
namespace transformer {
class LDBAProductBuilder {
   public:
    LDBAProductBuilder(const storm::automata::LimitDeterministicAutomaton& ldba, const std::vector<storm::storage::BitVector>& statesForAP)
        : ldba(ldba), statesForAP(statesForAP) {}

    template <typename Model>
    typename DAProduct<Model>::ptr build(const Model& originalMdp, const storm::storage::BitVector& statesOfInterest) const {
        return build<Model>(originalMdp.getTransitionMatrix(), originalMdp.getBackwardTransitions(), statesOfInterest);
    }

    template <typename Model>
    typename DAProduct<Model>::ptr build(const storm::storage::SparseMatrix<typename Model::ValueType>& transitionMatrix,
                                        const storm::storage::SparseMatrix<typename Model::ValueType>& backwardTransitions,
                                        const storm::storage::BitVector& statesOfInterest) const {
        storage::MaximalEndComponentDecomposition mecs(transitionMatrix, backwardTransitions);
        storm::storage::BitVector mecStates(transitionMatrix.getRowGroupCount());
        for (auto const& ec: mecs) {
            for (auto const& [state, _]: ec) {
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
        typename Product<Model>::ptr product = ProductBuilder<Model>::buildLimitDeterministicProduct(transitionMatrix, *this, statesOfInterest, mecStates);
        storm::automata::AcceptanceCondition::ptr prodAcceptance = ldba.getAcceptance()->lift(
            product->getProductModel().getNumberOfStates(), [&product](std::size_t prodState) { return product->getAutomatonState(prodState); });

        return typename DAProduct<Model>::ptr(new DAProduct<Model>(std::move(*product), prodAcceptance));
    }

    storm::storage::sparse::state_type getInitialState(storm::storage::sparse::state_type MdpState) const {
        return ldba.getInitialState();
    }

    const std::vector<storm::storage::sparse::state_type> getSuccessors(storm::storage::sparse::state_type automatonFrom, storm::storage::sparse::state_type modelTo) const {
        std::vector<storm::storage::sparse::state_type> succs;
        for (auto const& succ: ldba.getSuccessors(automatonFrom, getLabelForState(modelTo))) {
            succs.push_back(succ);
        }
        return succs;
    }

    const storm::storage::BitVector getAcceptingPart() const {
        return ldba.getAcceptingPart();
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

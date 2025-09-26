#pragma once

#include "storm/automata/Automaton.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/storage/BitVector.h"
#include "storm/transformer/DAProduct.h"
#include "storm/transformer/Product.h"
#include "storm/transformer/ProductBuilder.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"

#include <vector>

namespace storm::logic {
class PathFormula;
}
namespace storm {
namespace transformer {
class DAProductBuilder {
   public:
    DAProductBuilder(const storm::automata::DeterministicAutomaton& da, const std::vector<storm::storage::BitVector>& statesForAP)
        : da(da), statesForAP(statesForAP) {}

    template<typename ValueType>
    typename DAProduct<ValueType>::ptr build(storm::models::sparse::Model<ValueType> const& originalModel) const {
        // First build the product model which will be equipped with the initial states but no other label
        auto result = build(originalModel.getTransitionMatrix(), originalModel.getType(), originalModel.getInitialStates(), "init");
        liftOriginalModelInformationToProduct(originalModel, result);
        return result;
    }

    template<typename ValueType>
    typename DAProduct<ValueType>::ptr build(storm::storage::SparseMatrix<ValueType> const& originalMatrix, storm::models::ModelType const& modelType,
                                             storm::storage::BitVector const& statesOfInterest, std::string const& statesOfInterestLabel = "soi") const {
        typename Product<ValueType>::ptr product =
            ProductBuilder<ValueType>::buildProduct(originalMatrix, modelType, *this, statesOfInterest, statesOfInterestLabel);
        storm::automata::AcceptanceCondition::ptr prodAcceptance = da.getAcceptance()->lift(
            product->getProductModel().getNumberOfStates(), [&product](std::size_t prodState) { return product->getAutomatonState(prodState); });

        return typename DAProduct<ValueType>::ptr(new DAProduct<ValueType>(std::move(*product), prodAcceptance));
    }

    template<typename ValueType>
    static void liftOriginalModelInformationToProduct(storm::models::sparse::Model<ValueType> const& originalModel,
                                                      typename DAProduct<ValueType>::ptr product) {
        auto& productModel = product->getProductModel();
        // Then lift remaining model information to product model
        // State labeling
        for (auto const& label : originalModel.getStateLabeling().getLabels()) {
            if (label != "init") {  // init label already handled when building the product
                productModel.getStateLabeling().addLabel(label, product->liftFromModel(originalModel.getStateLabeling().getStates(label)));
            }
        }
        // Reward models
        for (auto const& [rewName, rewModel] : originalModel.getRewardModels()) {
            std::optional<std::vector<ValueType>> newStateRew, newStateActionRew;
            if (rewModel.hasStateRewards()) {
                newStateRew = product->liftFromModel(rewModel.getStateRewardVector());
            }
            if (rewModel.hasStateActionRewards()) {
                if (originalModel.isNondeterministicModel()) {
                    newStateActionRew =
                        product->liftChoiceVectorFromModel(rewModel.getStateActionRewardVector(), originalModel.getTransitionMatrix().getRowGroupIndices());
                } else {
                    newStateActionRew = product->liftFromModel(rewModel.getStateActionRewardVector());
                }
            }
            STORM_LOG_THROW(!rewModel.hasTransitionRewards(), storm::exceptions::NotSupportedException,
                            "Transition rewards are not supported in product models.");
            storm::models::sparse::StandardRewardModel<ValueType> productRewModel(std::move(newStateRew), std::move(newStateActionRew));
            productModel.getRewardModels().emplace(rewName, std::move(productRewModel));
        }
        // Continuous time information
        if (originalModel.isOfType(storm::models::ModelType::Ctmc)) {
            auto const& origCtmc = *originalModel.template as<storm::models::sparse::Ctmc<ValueType>>();
            auto& productCtmc = *productModel.template as<storm::models::sparse::Ctmc<ValueType>>();
            productCtmc.getExitRateVector() = product->liftFromModel(origCtmc.getExitRateVector());
        } else if (originalModel.isOfType(storm::models::ModelType::MarkovAutomaton)) {
            auto const& origMa = *originalModel.template as<storm::models::sparse::MarkovAutomaton<ValueType>>();
            auto& productMa = *productModel.template as<storm::models::sparse::MarkovAutomaton<ValueType>>();
            productMa.getMarkovianStates() = product->liftFromModel(origMa.getMarkovianStates());
            productMa.getExitRates() = product->liftChoiceVectorFromModel(origMa.getExitRates(), origMa.getTransitionMatrix().getRowGroupIndices());
        } else {
            STORM_LOG_ASSERT(originalModel.isOfType(storm::models::ModelType::Dtmc) || originalModel.isOfType(storm::models::ModelType::Mdp),
                             "Model type not supported for product construction");
        }
    }

    storm::storage::sparse::state_type getInitialState(storm::storage::sparse::state_type modelState) const {
        return da.getSuccessor(da.getInitialState(), getLabelForState(modelState));
    }

    storm::storage::sparse::state_type getSuccessor(storm::storage::sparse::state_type automatonFrom, storm::storage::sparse::state_type modelTo) const {
        return da.getSuccessor(automatonFrom, getLabelForState(modelTo));
    }

   private:
    const storm::automata::DeterministicAutomaton& da;
    const std::vector<storm::storage::BitVector>& statesForAP;

    storm::automata::APSet::alphabet_element getLabelForState(storm::storage::sparse::state_type s) const {
        storm::automata::APSet::alphabet_element label = da.getAPSet().elementAllFalse();
        for (unsigned int ap = 0; ap < da.getAPSet().size(); ap++) {
            if (statesForAP.at(ap).get(s)) {
                label = da.getAPSet().elementAddAP(label, ap);
            }
        }
        return label;
    }
};
}  // namespace transformer
}  // namespace storm

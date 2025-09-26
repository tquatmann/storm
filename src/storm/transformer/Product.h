#pragma once

#include <memory>
#include "storm/models/sparse/Model.h"

namespace storm {
namespace transformer {
template<typename ValueType>
class Product {
   public:
    typedef std::shared_ptr<Product<ValueType>> ptr;
    typedef storm::models::sparse::Model<ValueType> Model;
    typedef storm::storage::sparse::state_type state_type;
    typedef std::pair<state_type, state_type> product_state_type;
    typedef std::map<product_state_type, state_type> product_state_to_product_index_map;
    typedef std::vector<product_state_type> product_index_to_product_state_vector;

    Product(std::shared_ptr<Model> productModel, std::string&& productStateOfInterestLabel, product_state_to_product_index_map&& productStateToProductIndex,
            product_index_to_product_state_vector&& productIndexToProductState)
        : productModel(productModel),
          productStateOfInterestLabel(productStateOfInterestLabel),
          productStateToProductIndex(productStateToProductIndex),
          productIndexToProductState(productIndexToProductState) {}

    Product(Product<ValueType>&& product) = default;
    Product& operator=(Product<ValueType>&& product) = default;

    Model& getProductModel() {
        return *productModel;
    }

    std::shared_ptr<Model> getProductModelPtr() {
        return productModel;
    }

    state_type getModelState(state_type productStateIndex) const {
        return productIndexToProductState.at(productStateIndex).first;
    }

    state_type getAutomatonState(state_type productStateIndex) const {
        return productIndexToProductState.at(productStateIndex).second;
    }

    state_type getProductStateIndex(state_type modelState, state_type automatonState) const {
        return productStateToProductIndex.at(product_state_type(modelState, automatonState));
    }

    bool isValidProductState(state_type modelState, state_type automatonState) const {
        return (productStateToProductIndex.count(product_state_type(modelState, automatonState)) > 0);
    }

    storm::storage::BitVector liftFromAutomaton(const storm::storage::BitVector& vector) const {
        state_type n = productModel->getNumberOfStates();
        storm::storage::BitVector lifted(n, false);
        for (state_type s = 0; s < n; s++) {
            if (vector.get(getAutomatonState(s))) {
                lifted.set(s);
            }
        }
        return lifted;
    }

    storm::storage::BitVector liftFromModel(const storm::storage::BitVector& vector) const {
        state_type const n = productModel->getNumberOfStates();
        storm::storage::BitVector lifted(n, false);
        for (state_type s = 0; s < n; s++) {
            if (vector.get(getModelState(s))) {
                lifted.set(s);
            }
        }
        return lifted;
    }

    template<typename T>
    std::vector<T> liftFromModel(const std::vector<T>& vector) const {
        state_type const n = productModel->getNumberOfStates();
        std::vector<T> lifted;
        lifted.reserve(n);
        for (state_type s = 0; s < n; s++) {
            lifted.push_back(vector[getModelState(s)]);
        }
        return lifted;
    }

    template<typename T>
    std::vector<T> liftChoiceVectorFromModel(const std::vector<T>& vector, std::vector<uint64_t> const& modelRowGroupIndices) const {
        state_type const n = productModel->getNumberOfStates();
        std::vector<T> lifted;
        lifted.reserve(productModel->getNumberOfChoices());
        for (state_type s = 0; s < n; s++) {
            auto const modelState = getModelState(s);
            for (auto choice = modelRowGroupIndices[modelState]; choice < modelRowGroupIndices[modelState + 1]; ++choice) {
                lifted.push_back(vector[choice]);
            }
        }
        return lifted;
    }

    template<typename T>
    std::vector<T> projectToOriginalModel(const Model& originalModel, const std::vector<T>& prodValues) {
        return projectToOriginalModel(originalModel.getNumberOfStates(), prodValues);
    }

    template<typename T>
    std::vector<T> projectToOriginalModel(std::size_t numberOfStates, const std::vector<T>& prodValues) {
        std::vector<T> origValues(numberOfStates);
        for (state_type productState : productModel->getStateLabeling().getStates(productStateOfInterestLabel)) {
            state_type originalState = getModelState(productState);
            origValues.at(originalState) = prodValues.at(productState);
        }
        return origValues;
    }

    const storm::storage::BitVector& getStatesOfInterest() const {
        return productModel->getStates(productStateOfInterestLabel);
    }

    void printMapping(std::ostream& out) const {
        out << "Mapping index -> product state\n";
        for (std::size_t i = 0; i < productIndexToProductState.size(); i++) {
            out << " " << i << ": " << productIndexToProductState.at(i).first << "," << productIndexToProductState.at(i).second << "\n";
        }
    }

   private:
    std::shared_ptr<Model> productModel;
    std::string productStateOfInterestLabel;
    product_state_to_product_index_map productStateToProductIndex;
    product_index_to_product_state_vector productIndexToProductState;
};
}  // namespace transformer
}  // namespace storm

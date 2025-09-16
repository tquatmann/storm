#pragma once

#include <memory>
#include "storm/automata/AcceptanceCondition.h"
#include "storm/transformer/Product.h"

namespace storm {
namespace transformer {

template<typename Model>
class DAMultiProduct : public Product<Model> {
   public:
    typedef std::shared_ptr<DAMultiProduct<Model>> ptr;

    DAMultiProduct(Product<Model>&& product, std::vector<storm::automata::AcceptanceCondition::ptr> acceptances)
        : Product<Model>(std::move(product)), acceptances(acceptances) {
        // Intentionally left blank
    }

    std::vector<storm::automata::AcceptanceCondition::ptr> getAcceptances() {
        return acceptances;
    }

   private:
    std::vector<storm::automata::AcceptanceCondition::ptr> acceptances;
};
}  // namespace transformer
}  // namespace storm

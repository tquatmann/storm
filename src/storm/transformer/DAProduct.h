#pragma once

#include <memory>
#include "storm/automata/AcceptanceCondition.h"
#include "storm/transformer/Product.h"

namespace storm {
namespace transformer {

template<typename ValueType>
class DAProduct : public Product<ValueType> {
   public:
    typedef std::shared_ptr<DAProduct<ValueType>> ptr;

    DAProduct(Product<ValueType>&& product, storm::automata::AcceptanceCondition::ptr acceptance)
        : Product<ValueType>(std::move(product)), acceptance(acceptance) {
        // Intentionally left blank
    }

    storm::automata::AcceptanceCondition::ptr getAcceptance() {
        return acceptance;
    }

   private:
    storm::automata::AcceptanceCondition::ptr acceptance;
};
}  // namespace transformer
}  // namespace storm

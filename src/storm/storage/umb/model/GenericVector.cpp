#include "storm/storage/umb/model/GenericVector.h"

#include "storm/storage/umb/model/ValueEncoding.h"
namespace storm::umb {

void GenericVector::unset() {
    data = std::monostate();
}

bool GenericVector::hasValue() const {
    return !std::holds_alternative<std::monostate>(data);
}

}  // namespace storm::umb
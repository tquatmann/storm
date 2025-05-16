#pragma once

#include <variant>

#include "storm/storage/umb/model/StorageType.h"
#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {
class UmbModelBase {
   public:
    template<StorageType Storage>
    using UmbModelPtr = std::unique_ptr<UmbModel<Storage>>;

    template<StorageType Storage>
    UmbModelBase(UmbModelPtr<Storage>&& model) : data(std::move(model)) {}

    bool isStorageType(StorageType storageType) const;
    ModelIndex const& getIndex() const;
    bool validate(bool silent = false) const;

    template<StorageType Storage>
    UmbModel<Storage> const& as() const;

    template<StorageType Storage>
        requires(Storage == StorageType::Memory)
    UmbModel<Storage>& as();

   private:
    // Note: we want to use boost::pfr for UmbModels, which means that proper polymorphism is not possible
    // This class is a poor man's variant implementation instead of inheritance
    std::variant<UmbModelPtr<StorageType::Disk>, UmbModelPtr<StorageType::Memory>> data;
};

}  // namespace storm::umb

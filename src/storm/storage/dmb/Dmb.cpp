#include "storm/storage/dmb/Dmb.h"

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/storage/dmb/export/DmbExport.h"
#include "storm/storage/dmb/export/SparseModelToDmb.h"
#include "storm/storage/dmb/import/DmbImport.h"
#include "storm/storage/dmb/import/SparseModelFromDmb.h"

namespace storm::dmb {

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModelFromDmb(std::filesystem::path const& dmbLocation,
                                                                                            ImportOptions const& options) {
    auto dmb = storm::dmb::fromDisk(dmbLocation, options);
    return storm::dmb::sparseModelFromDmb<ValueType, RewardModelType>(*dmb);
}

template<typename ValueType, typename RewardModelType>
void exportModelToDmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model, std::filesystem::path const& targetLocation,
                      ExportOptions const& options) {
    auto dmb = storm::dmb::sparseModelToDmb(model);
    storm::dmb::toDisk(*dmb, targetLocation, options);
}

template std::shared_ptr<storm::models::sparse::Model<double>> parseModelFromDmb(std::filesystem::path const&, ImportOptions const&);
template std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> parseModelFromDmb(std::filesystem::path const&, ImportOptions const&);

template void exportModelToDmb(storm::models::sparse::Model<double> const&, std::filesystem::path const&, ExportOptions const&);
template void exportModelToDmb(storm::models::sparse::Model<storm::RationalNumber> const&, std::filesystem::path const&, ExportOptions const&);

}  // namespace storm::dmb

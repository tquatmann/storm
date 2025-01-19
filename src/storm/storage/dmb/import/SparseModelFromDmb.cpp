#include "storm/storage/dmb/import/SparseModelFromDmb.h"

#include "storm/storage/dmb/model/DmbModel.h"

#include "storm/models/ModelType.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/Smg.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/sparse/ModelComponents.h"
#include "storm/utility/builder.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"

namespace storm::dmb {

namespace detail {

template<typename ValueType, StorageType Storage>
storm::storage::SparseMatrix<ValueType> constructTransitionMatrix(storm::dmb::DmbModel<Storage> const& dmbModel) {
    auto const& index = dmbModel.getIndex();
    bool const hasRowGroups = dmbModel.states.stateToChoice.has_value();
    storm::storage::SparseMatrixBuilder<ValueType> builder(index.numChoices, index.numStates, index.numBranches, true, hasRowGroups,
                                                           hasRowGroups ? index.numStates : 0u);
    for (uint64_t stateIndex{0}, choiceIndex{0}, branchIndex{0}; stateIndex < index.numStates; ++stateIndex) {
        if (hasRowGroups) {
            builder.newRowGroup(choiceIndex);
        }
        for (auto const choiceEnd = hasRowGroups ? dmbModel.states.stateToChoice.value()[stateIndex + 1] : stateIndex + 1; choiceIndex < choiceEnd;
             ++choiceIndex) {
            STORM_LOG_ASSERT(choiceIndex < index.numChoices, "Choice index out of bounds.");
            for (auto const branchEnd = dmbModel.choices.choiceToBranch.value()[choiceIndex + 1]; branchIndex < branchEnd; ++branchIndex) {
                STORM_LOG_ASSERT(branchIndex < index.numBranches, "branch index out of bounds.");
                builder.addNextValue(choiceIndex, dmbModel.branches.branchToTarget.value()[branchIndex],
                                     dmbModel.branches.branchToValue.template at<ValueType>(branchIndex));
            }
        };
    }
    return builder.build();
}

template<typename ValueType, typename RewardModelType, StorageType Storage>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> constructSparseModel(storm::dmb::DmbModel<Storage> const& dmbModel) {
    storm::storage::sparse::ModelComponents<ValueType, RewardModelType> components(constructTransitionMatrix<ValueType, Storage>(dmbModel));
    return storm::utility::builder::buildModelFromComponents(deriveModelType(dmbModel.getIndex()), std::move(components));
}

}  // namespace detail

storm::models::ModelType deriveModelType(storm::dmb::ModelIndex const& index) {
    STORM_LOG_THROW(index.branchValues != storm::dmb::ModelIndex::BranchValues::None, storm::exceptions::NotSupportedException,
                    "Models without branch values are not supported.");
    switch (index.time) {
        using enum storm::models::ModelType;
        using enum storm::dmb::ModelIndex::Time;
        case Discrete:
            switch (index.players) {
                case 0:
                    return Dtmc;
                case 1:
                    return Mdp;
                default:
                    return Smg;
            }
        case Stochastic:
            STORM_LOG_THROW(index.players == 0, storm::exceptions::NotSupportedException, "Stochastic time models with multiple players are not supported.");
            return Ctmc;
        case UrgentStochastic:
            STORM_LOG_THROW(index.players == 1, storm::exceptions::NotSupportedException,
                            "Urgent stochastic time models with multiple or no players are not supported.");
            return MarkovAutomaton;
    }
}

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> sparseModelFromDmb(storm::dmb::DmbModelBase const& dmbModel) {
    if (dmbModel.isStorageType(StorageType::Disk)) {
        return detail::constructSparseModel<ValueType, RewardModelType, StorageType::Disk>(dmbModel.as<StorageType::Disk>());
    } else if (dmbModel.isStorageType(StorageType::Memory)) {
        return detail::constructSparseModel<ValueType, RewardModelType, StorageType::Memory>(dmbModel.as<StorageType::Memory>());
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Storage type not expected.");
    }
}
template std::shared_ptr<storm::models::sparse::Model<double>> sparseModelFromDmb<double>(storm::dmb::DmbModelBase const& dmbModel);
template std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> sparseModelFromDmb<storm::RationalNumber>(
    storm::dmb::DmbModelBase const& dmbModel);
}  // namespace storm::dmb
#include "storm/storage/umb/import/SparseModelFromUmb.h"

#include "storm/storage/umb/model/UmbModel.h"

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

namespace storm::umb {

namespace detail {

template<typename ValueType, StorageType Storage>
storm::storage::SparseMatrix<ValueType> constructTransitionMatrix(storm::umb::UmbModel<Storage> const& umbModel) {
    auto const& ts = umbModel.getIndex().transitionSystem;
    bool const hasRowGroups = umbModel.states.stateToChoice.has_value();
    storm::storage::SparseMatrixBuilder<ValueType> builder(ts.numChoices, ts.numStates, ts.numBranches, true, hasRowGroups, hasRowGroups ? ts.numStates : 0u);
    for (uint64_t stateIndex{0}, choiceIndex{0}, branchIndex{0}; stateIndex < ts.numStates; ++stateIndex) {
        if (hasRowGroups) {
            builder.newRowGroup(choiceIndex);
        }
        for (auto const choiceEnd = hasRowGroups ? umbModel.states.stateToChoice.value()[stateIndex + 1] : stateIndex + 1; choiceIndex < choiceEnd;
             ++choiceIndex) {
            STORM_LOG_ASSERT(choiceIndex < ts.numChoices, "Choice index out of bounds.");
            for (auto const branchEnd = umbModel.choices.choiceToBranch.value()[choiceIndex + 1]; branchIndex < branchEnd; ++branchIndex) {
                STORM_LOG_ASSERT(branchIndex < ts.numBranches, "branch index out of bounds.");
                builder.addNextValue(choiceIndex, umbModel.branches.branchToTarget.value()[branchIndex],
                                     umbModel.branches.branchToValue.template at<ValueType>(branchIndex));
            }
        };
    }
    return builder.build();
}

template<StorageType Storage>
storm::storage::BitVector createBitVector(storm::umb::VectorType<bool, Storage> const& umbBitVector, uint64_t size) {
    if constexpr (Storage == StorageType::Memory) {
        STORM_LOG_ASSERT(umbBitVector.size() >= size, "Bit vector has unexpected size: " << umbBitVector.size() << " < " << size);
        storm::storage::BitVector result = umbBitVector;
        result.resize(size);
        return result;
    } else {
        return umbBitVector.getAsBitVector(size);
    }
}

template<StorageType Storage>
storm::models::sparse::StateLabeling constructStateLabelling(storm::umb::UmbModel<Storage> const& umbModel) {
    auto const& ts = umbModel.getIndex().transitionSystem;
    storm::models::sparse::StateLabeling stateLabelling(ts.numStates);
    if (umbModel.states.initialStates) {
        stateLabelling.addLabel("init", createBitVector<Storage>(umbModel.states.initialStates.value(), ts.numStates));
    }
    // todo: more labels
    return stateLabelling;
}

template<typename ValueType, typename RewardModelType, StorageType Storage>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> constructSparseModel(storm::umb::UmbModel<Storage> const& umbModel) {
    storm::storage::sparse::ModelComponents<ValueType, RewardModelType> components(constructTransitionMatrix<ValueType, Storage>(umbModel),
                                                                                   constructStateLabelling(umbModel));
    return storm::utility::builder::buildModelFromComponents(deriveModelType(umbModel.getIndex()), std::move(components));
}

}  // namespace detail

storm::models::ModelType deriveModelType(storm::umb::ModelIndex const& index) {
    auto const& ts = index.transitionSystem;

    STORM_LOG_THROW(ts.branchValues != storm::umb::ModelIndex::TransitionSystem::BranchValues::None, storm::exceptions::NotSupportedException,
                    "Models without branch values are not supported.");
    switch (ts.time) {
        using enum storm::models::ModelType;
        using enum storm::umb::ModelIndex::TransitionSystem::Time;
        case Discrete:
            switch (ts.players) {
                case 0:
                    return Dtmc;
                case 1:
                    return Mdp;
                default:
                    return Smg;
            }
        case Stochastic:
            STORM_LOG_THROW(ts.players == 0, storm::exceptions::NotSupportedException, "Stochastic time models with multiple players are not supported.");
            return Ctmc;
        case UrgentStochastic:
            STORM_LOG_THROW(ts.players == 1, storm::exceptions::NotSupportedException,
                            "Urgent stochastic time models with multiple or no players are not supported.");
            return MarkovAutomaton;
    }
}

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> sparseModelFromUmb(storm::umb::UmbModelBase const& umbModel) {
    if (umbModel.isStorageType(StorageType::Disk)) {
        return detail::constructSparseModel<ValueType, RewardModelType, StorageType::Disk>(umbModel.as<StorageType::Disk>());
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        return detail::constructSparseModel<ValueType, RewardModelType, StorageType::Memory>(umbModel.as<StorageType::Memory>());
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Storage type not expected.");
    }
}
template std::shared_ptr<storm::models::sparse::Model<double>> sparseModelFromUmb<double>(storm::umb::UmbModelBase const& umbModel);
template std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> sparseModelFromUmb<storm::RationalNumber>(
    storm::umb::UmbModelBase const& umbModel);
}  // namespace storm::umb
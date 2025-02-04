#include "storm/storage/dmb/export/SparseModelToDmb.h"

#include "storm/storage/SparseMatrix.h"
#include "storm/storage/dmb/model/DmbModel.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/models/ModelType.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/Smg.h"
#include "storm/utility/macros.h"

namespace storm::dmb {

namespace detail {

template<typename ValueType>
void transitionMatrixToDmb(storm::storage::SparseMatrix<ValueType> const& matrix, storm::dmb::DmbModel<StorageType::Memory>& dmb) {
    if (!matrix.hasTrivialRowGrouping()) {
        dmb.states.stateToChoice = matrix.getRowGroupIndices();
    }
    dmb.choices.choiceToBranch = matrix.getRowIndices();
    std::vector<uint64_t> branchToTarget;
    branchToTarget.reserve(matrix.getEntryCount());
    std::vector<ValueType> branchToValue;
    branchToValue.reserve(matrix.getEntryCount());
    for (auto const& entry : matrix) {
        branchToTarget.push_back(entry.getColumn());
        branchToValue.push_back(entry.getValue());
    }
    dmb.branches.branchToTarget = std::move(branchToTarget);
    dmb.branches.branchToValue.template set<ValueType>(std::move(branchToValue));
}

template<typename ValueType, typename RewardModelType>
void setIndexInformation(storm::models::sparse::Model<ValueType, RewardModelType> const& model, storm::dmb::ModelIndex& index) {
    index.creation.setDateToNow();
    index.creation.tool = "storm";
    // TODO: it's apparently difficult to get the version of the tool because the storm library is not linked against storm-version-info

    index.branchValues = storm::dmb::ModelIndex::BranchValues::Number;
    if (std::is_same_v<ValueType, double>) {
        index.branchValueType = storm::dmb::ModelIndex::BranchValueType::Double;
    } else if (std::is_same_v<ValueType, storm::RationalNumber>) {
        index.branchValueType = storm::dmb::ModelIndex::BranchValueType::Rational;
    } else if (std::is_same_v<ValueType, storm::Interval>) {
        index.branchValueType = storm::dmb::ModelIndex::BranchValueType::Double;
        index.branchValues = storm::dmb::ModelIndex::BranchValues::Interval;
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected ValueType.");
    }

    index.numStates = model.getNumberOfStates();
    index.numInitialStates = model.getInitialStates().getNumberOfSetBits();
    index.numChoices = model.getNumberOfChoices();
    index.numBranches = model.getNumberOfTransitions();

    switch (model.getType()) {
        using enum storm::models::ModelType;
        using enum storm::dmb::ModelIndex::Time;
        case Dtmc:
            index.time = Discrete;
            index.players = 0;
            break;
        case Ctmc:
            index.time = Stochastic;
            index.players = 0;
            break;
        case Mdp:
        case Pomdp:
            index.time = Discrete;
            index.players = 1;
            break;
        case MarkovAutomaton:
            index.time = UrgentStochastic;
            index.players = 1;
            break;
        case Smg:
            index.time = Discrete;
            index.players = dynamic_cast<storm::models::sparse::Smg<ValueType, RewardModelType> const&>(model).getNumberOfPlayers();
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected model type.");
    }
    index.time = storm::dmb::ModelIndex::Time::Discrete;
    index.branchValues = storm::dmb::ModelIndex::BranchValues::Number;
}
}  // namespace detail

template<typename ValueType, typename RewardModelType>
std::unique_ptr<storm::dmb::DmbModelBase> sparseModelToDmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model) {
    DmbModel<StorageType::Memory> dmbModel;
    detail::setIndexInformation(model, dmbModel.index);
    if (model.isOfType(storm::models::ModelType::Ctmc)) {
        detail::transitionMatrixToDmb(dynamic_cast<storm::models::sparse::Ctmc<ValueType> const&>(model).computeProbabilityMatrix(), dmbModel);
    } else {
        detail::transitionMatrixToDmb(model.getTransitionMatrix(), dmbModel);
    }
    return std::make_unique<DmbModel<StorageType::Memory>>(std::move(dmbModel));
}
template std::unique_ptr<storm::dmb::DmbModelBase> sparseModelToDmb<double>(storm::models::sparse::Model<double> const& model);
template std::unique_ptr<storm::dmb::DmbModelBase> sparseModelToDmb<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model);
}  // namespace storm::dmb
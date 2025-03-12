#include "storm/storage/umb/export/SparseModelToUmb.h"

#include "storm/storage/SparseMatrix.h"
#include "storm/storage/umb/model/UmbModel.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/WrongFormatException.h"
#include "storm/models/ModelType.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/Smg.h"
#include "storm/utility/macros.h"

namespace storm::umb {

namespace detail {

bool isValidAnnotationId(std::string const& id) {
    return !id.empty() && std::all_of(id.begin(), id.end(), [](char c) { return std::isalnum(c) || c == '_' || c == '-'; });
}

template<typename ValueType>
void transitionMatrixToUmb(storm::storage::SparseMatrix<ValueType> const& matrix, storm::umb::UmbModel<StorageType::Memory>& umb) {
    if (!matrix.hasTrivialRowGrouping()) {
        umb.states.stateToChoice = matrix.getRowGroupIndices();
    }
    umb.choices.choiceToBranch = matrix.getRowIndices();
    std::vector<uint64_t> branchToTarget;
    branchToTarget.reserve(matrix.getEntryCount());
    std::vector<ValueType> branchToValue;
    branchToValue.reserve(matrix.getEntryCount());
    for (auto const& entry : matrix) {
        branchToTarget.push_back(entry.getColumn());
        branchToValue.push_back(entry.getValue());
    }
    umb.branches.branchToTarget = std::move(branchToTarget);
    if constexpr (std::is_same_v<ValueType, double>) {
        umb.branches.branchValues.template set<ValueType>(std::move(branchToValue));
    } else if constexpr (std::is_same_v<ValueType, storm::RationalNumber>) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Rationals not yet supported.");
    } else {
        static_assert(false, "Unhandled value type.");
    }
}

void stateLabelingToUmb(storm::models::sparse::StateLabeling const& labeling, storm::umb::UmbModel<StorageType::Memory>& umb) {
    for (auto const& labelName : labeling.getLabels()) {
        if (labelName == "init") {
            umb.states.initialStates = labeling.getStates(labelName);
        } else {
            // Get a unique, valid id
            auto const labelIdBase = isValidAnnotationId(labelName) ? labelName : "label";
            std::string labelId = labelIdBase;
            for (uint64_t i = 0; umb.annotations.contains(labelId); ++i) {
                labelId = labelIdBase + "_" + std::to_string(i);
            }

            // Add the annotation data
            umb.annotations[labelId].values.template set<bool>(labeling.getStates(labelName));

            // Add index data
            STORM_LOG_ASSERT(!umb.index.annotations.contains(labelId), "Annotation id " << labelId << " already exists in index but not as file.");
            auto& annotationIndex = umb.index.annotations[labelId];
            annotationIndex.name = labelName;
            annotationIndex.appliesTo = storm::umb::ModelIndex::Annotation::AppliesTo::States;
            annotationIndex.type = storm::umb::ModelIndex::Annotation::Type::Bool;
        }
    }
}

template<typename ValueType, typename RewardModelType>
void setIndexInformation(storm::models::sparse::Model<ValueType, RewardModelType> const& model, storm::umb::ModelIndex& index) {
    index.fileData.emplace();
    index.fileData->setDateToNow();
    index.fileData->tool = "storm";
    // TODO: it's apparently difficult to get the version of the tool because the storm library is not linked against storm-version-info

    auto& ts = index.transitionSystem;
    {
        using enum storm::umb::ModelIndex::TransitionSystem::BranchValues;
        using enum storm::umb::ModelIndex::TransitionSystem::BranchValueType;
        ts.branchValues = Number;
        if (std::is_same_v<ValueType, double>) {
            ts.branchValueType = Double;
        } else if (std::is_same_v<ValueType, storm::RationalNumber>) {
            ts.branchValueType = Rational;
        } else if (std::is_same_v<ValueType, storm::Interval>) {
            ts.branchValueType = Double;
            ts.branchValues = Interval;
        } else {
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected ValueType.");
        }
    }

    ts.numStates = model.getNumberOfStates();
    ts.numInitialStates = model.getInitialStates().getNumberOfSetBits();
    ts.numChoices = model.getNumberOfChoices();
    ts.numBranches = model.getNumberOfTransitions();

    switch (model.getType()) {
        using enum storm::models::ModelType;
        using enum storm::umb::ModelIndex::TransitionSystem::Time;
        case Dtmc:
            ts.time = Discrete;
            ts.players = 0;
            break;
        case Ctmc:
            ts.time = Stochastic;
            ts.players = 0;
            break;
        case Mdp:
        case Pomdp:
            ts.time = Discrete;
            ts.players = 1;
            break;
        case MarkovAutomaton:
            ts.time = UrgentStochastic;
            ts.players = 1;
            break;
        case Smg:
            ts.time = Discrete;
            ts.players = dynamic_cast<storm::models::sparse::Smg<ValueType, RewardModelType> const&>(model).getNumberOfPlayers();
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected model type.");
    }
}
}  // namespace detail

template<typename ValueType, typename RewardModelType>
std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb(storm::models::sparse::Model<ValueType, RewardModelType> const& model) {
    UmbModel<StorageType::Memory> umbModel;
    detail::setIndexInformation(model, umbModel.index);
    detail::stateLabelingToUmb(model.getStateLabeling(), umbModel);
    if (model.isOfType(storm::models::ModelType::Ctmc)) {
        detail::transitionMatrixToUmb(dynamic_cast<storm::models::sparse::Ctmc<ValueType> const&>(model).computeProbabilityMatrix(), umbModel);
    } else {
        detail::transitionMatrixToUmb(model.getTransitionMatrix(), umbModel);
    }
    return std::make_unique<UmbModel<StorageType::Memory>>(std::move(umbModel));
}

template std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb<double>(storm::models::sparse::Model<double> const& model);
template std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model);
}  // namespace storm::umb
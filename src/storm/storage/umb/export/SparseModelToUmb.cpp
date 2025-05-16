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
#include "storm/utility/vector.h"

namespace storm::umb {

namespace detail {

bool isValidAnnotationId(std::string const& id) {
    return !id.empty() && std::all_of(id.begin(), id.end(), [](char c) { return std::isalnum(c) || c == '_' || c == '-'; });
}

template<typename ValueType, typename TargetValueType>
void transitionMatrixToUmb(storm::storage::SparseMatrix<ValueType> const& matrix, storm::umb::UmbModel<StorageType::Memory>& umb, bool normalize) {
    if (!matrix.hasTrivialRowGrouping()) {
        umb.states.stateToChoice = matrix.getRowGroupIndices();
    }
    umb.choices.choiceToBranch = matrix.getRowIndices();
    umb.branches.branchToTarget.emplace().reserve(matrix.getEntryCount());
    std::vector<TargetValueType> branchProbabilities;
    branchProbabilities.reserve(matrix.getEntryCount());
    for (uint64_t rowIndex = 0; rowIndex < matrix.getRowCount(); ++rowIndex) {
        auto const& row = matrix.getRow(rowIndex);
        for (auto const& entry : row) {
            umb.branches.branchToTarget->push_back(entry.getColumn());
            branchProbabilities.push_back(storm::utility::convertNumber<TargetValueType>(entry.getValue()));
        }
        if (normalize) {
            auto rowProbs = std::span<TargetValueType>(branchProbabilities.end() - row.getNumberOfEntries(), branchProbabilities.end());
            TargetValueType const rowSum = std::accumulate(rowProbs.begin(), rowProbs.end(), storm::utility::zero<TargetValueType>());
            if (!storm::utility::isOne(rowSum)) {
                std::for_each(rowProbs.begin(), rowProbs.end(), [&rowSum](TargetValueType& entry) { entry /= rowSum; });
            }
        }
    }
    umb.branches.branchProbabilities.template set<TargetValueType>(std::move(branchProbabilities));
}

void labelingToUmb(storm::models::sparse::ItemLabeling const& labeling, storm::umb::UmbModel<StorageType::Memory>& umb) {
    for (auto const labelName : labeling.getLabels()) {
        if (labeling.isStateLabeling() && labelName == "init") {
            continue;  // skip initial state labeling.
        }
        auto const name = umb.index.annotations.findAtomicPropositionName(labelName);
        STORM_LOG_ASSERT(!name.empty(), "Label '" << labelName << "' not found in the model index.");
        auto& annotation = umb.atomicPropositions[name];
        if (labeling.isStateLabeling()) {
            annotation.forStates.emplace().values = labeling.asStateLabeling().getStates(labelName);
        } else {
            STORM_LOG_ASSERT(labeling.isChoiceLabeling(), "Unexpected item label.");
            annotation.forChoices.emplace().values = labeling.asChoiceLabeling().getChoices(labelName);
        }
    }
}

template<typename ValueType, typename TargetValueType>
void rewardToUmb(std::string const& identifier, storm::models::sparse::StandardRewardModel<ValueType> const& rewardModel,
                 storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::umb::UmbModel<StorageType::Memory>& umb) {
    auto const rewardName = umb.index.annotations.findRewardName(identifier);
    STORM_LOG_ASSERT(!rewardName.empty(), "Reward '" << identifier << "' not found in the model index.");
    // auto const& rewardIndex = umb.index.annotations.rewards.at(rewardName);
    STORM_LOG_ASSERT(!umb.rewards.contains(rewardName), "Reward '" << identifier << "' already exists in the umb model.");
    auto& rewardAnnotation = umb.rewards[rewardName];
    if (rewardModel.hasStateRewards()) {
        rewardAnnotation.forStates.emplace().values.template set<TargetValueType>(
            storm::utility::vector::convertNumericVector<TargetValueType>(rewardModel.getStateRewardVector()));
    }
    if (rewardModel.hasStateActionRewards()) {
        rewardAnnotation.forChoices.emplace().values.template set<TargetValueType>(
            storm::utility::vector::convertNumericVector<TargetValueType>(rewardModel.getStateRewardVector()));
    }
    if (rewardModel.hasTransitionRewards()) {
        std::vector<TargetValueType> branchRewards;
        branchRewards.reserve(transitionMatrix.getEntryCount());
        STORM_LOG_ASSERT(transitionMatrix.getRowCount() == rewardModel.getTransitionRewardMatrix().getRowCount(),
                         "The number of rows in the transition matrix and the reward model do not match.");
        for (uint64_t rowIndex = 0; rowIndex < transitionMatrix.getRowCount(); ++rowIndex) {
            auto const& transitionRow = transitionMatrix.getRow(rowIndex);
            auto const& rewardRow = rewardModel.getTransitionRewardMatrix().getRow(rowIndex);
            auto rewIt = rewardRow.begin();
            // Match transition branch entries with entries in the transition reward matrix (which might not have the same entries at the same columns)
            for (auto const& entry : transitionRow) {
                while (rewIt != rewardRow.end() && rewIt->getColumn() < entry.getColumn()) {
                    ++rewIt;
                }
                if (rewIt == rewardRow.end() || rewIt->getColumn() > entry.getColumn()) {
                    branchRewards.push_back(storm::utility::zero<TargetValueType>());
                } else {
                    STORM_LOG_ASSERT(rewIt->getColumn() == entry.getColumn(), "Unexpected column in reward model.");
                    branchRewards.push_back(storm::utility::convertNumber<TargetValueType>(rewIt->getValue()));
                }
            }
        }
        rewardAnnotation.forBranches.emplace().values.template set<TargetValueType>(std::move(branchRewards));
    }
}

template<typename ValueType>
void setIndexInformation(storm::models::sparse::Model<ValueType> const& model, storm::umb::ModelIndex& index, ExportOptions const& options) {
    // No model (meta-)data to set at this point.
    // file-data:
    index.fileData.emplace();
    index.fileData->setCreationDateToNow();
    index.fileData->tool = "storm";
    // TODO: it's apparently difficult to get the version of the tool because the storm library is not linked against storm-version-info

    // transition-system:
    auto& ts = index.transitionSystem;
    switch (model.getType()) {
        using enum storm::models::ModelType;
        using enum storm::umb::ModelIndex::TransitionSystem::Time;
        case Dtmc:
            ts.time = Discrete;
            ts.numPlayers = 0;
            break;
        case Ctmc:
            ts.time = Stochastic;
            ts.numPlayers = 0;
            break;
        case Mdp:
        case Pomdp:
            ts.time = Discrete;
            ts.numPlayers = 1;
            break;
        case MarkovAutomaton:
            ts.time = UrgentStochastic;
            ts.numPlayers = 1;
            break;
        case Smg:
            ts.time = Discrete;
            ts.numPlayers = dynamic_cast<storm::models::sparse::Smg<ValueType> const&>(model).getNumberOfPlayers();
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected model type.");
    }
    ts.numStates = model.getNumberOfStates();
    ts.numInitialStates = model.getInitialStates().getNumberOfSetBits();
    ts.numChoices = model.getNumberOfChoices();
    ts.numActions = model.hasChoiceOrigins() ? model.getChoiceOrigins()->getNumberOfIdentifiers() : 1;
    ts.numBranches = model.getNumberOfTransitions();
    auto targetType = options.valueType;
    if (targetType == ExportOptions::ValueType::Default) {
        if (std::is_same_v<ValueType, double>) {
            targetType = ExportOptions::ValueType::Double;
        } else if (std::is_same_v<ValueType, storm::RationalNumber>) {
            targetType = ExportOptions::ValueType::Rational;
        } else if (std::is_same_v<ValueType, storm::Interval>) {
            targetType = ExportOptions::ValueType::DoubleInterval;
        } else {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected value type.");
        }
    }
    using Prob = storm::umb::ModelIndex::TransitionSystem::BranchProbabilityType;
    using Rew = storm::umb::ModelIndex::Annotations::Reward::Type;
    Rew rewardType;
    switch (targetType) {
        case ExportOptions::ValueType::Double:
            ts.branchProbabilityType = Prob::Double;
            rewardType = Rew::Double;
            break;
        case ExportOptions::ValueType::Rational:
            ts.branchProbabilityType = Prob::Rational;
            rewardType = Rew::Rational;
            break;
        case ExportOptions::ValueType::DoubleInterval:
            ts.branchProbabilityType = Prob::DoubleInterval;
            rewardType = Rew::DoubleInterval;
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected value type.");
    }

    // annotations:
    // rewards:
    auto& rewards = index.annotations.rewards;
    for (auto const& [rewardModelName, rewardModel] : model.getRewardModels()) {
        auto const [name, alias] = index.annotations.getAllowedNameAndAlias(rewardModelName);
        STORM_LOG_THROW(!rewards.contains(name), storm::exceptions::WrongFormatException, "Reward id '" << name << "' already exists.");
        auto& rewardIndex = rewards[name];
        rewardIndex.alias = alias;
        rewardIndex.type = rewardType;
        rewardIndex.appliesTo.emplace();
        using enum storm::umb::ModelIndex::Annotations::AppliesTo;
        if (rewardModel.hasStateRewards()) {
            rewardIndex.appliesTo->push_back(States);
        }
        if (rewardModel.hasStateActionRewards()) {
            rewardIndex.appliesTo->push_back(Choices);
        }
        if (rewardModel.hasTransitionRewards()) {
            rewardIndex.appliesTo->push_back(Branches);
        }
    }

    // aps:
    auto& aps = index.annotations.atomicPropositions;
    for (auto const& label : model.getStateLabeling().getLabels()) {
        if (label == "init") {
            continue;
        }
        auto const [name, alias] = index.annotations.getAllowedNameAndAlias(label);
        STORM_LOG_THROW(!aps.contains(name), storm::exceptions::WrongFormatException, "AP for state label '" << label << "' already exists.");
        auto& apIndex = aps[name];
        apIndex.alias = alias;
        apIndex.appliesTo.emplace();
        apIndex.appliesTo->push_back(storm::umb::ModelIndex::Annotations::AppliesTo::States);
    }
    if (model.hasChoiceLabeling()) {
        for (auto const& label : model.getChoiceLabeling().getLabels()) {
            auto const [name, alias] = index.annotations.getAllowedNameAndAlias(label);
            // Note: An AP might exist as both, choice and state label
            bool const existsAlready = aps.contains(name);
            auto& apIndex = aps[name];
            if (existsAlready) {
                // Assert that the labelling only exists as state label
                STORM_LOG_THROW(apIndex.appliesTo->size() == 1 && apIndex.appliesTo->front() == storm::umb::ModelIndex::Annotations::AppliesTo::States,
                                storm::exceptions::WrongFormatException, "AP for choice label '" << label << "' already exists.");
                // Assert that the alias match
                STORM_LOG_THROW(
                    apIndex.alias == alias, storm::exceptions::WrongFormatException,
                    "AP '" << name << "' exists under two different alias: '" << apIndex.alias.value_or("<none>") << " and " << alias.value_or("<none>"));
            } else {
                apIndex.alias = alias;
                apIndex.appliesTo.emplace();
            }
            apIndex.appliesTo->push_back(storm::umb::ModelIndex::Annotations::AppliesTo::Choices);
        }
    }
}

template<typename ValueType, typename TargetValueType>
void sparseModelToUmb(storm::models::sparse::Model<ValueType> const& model, UmbModel<StorageType::Memory>& umbModel, ExportOptions const& /* options */) {
    umbModel.states.initialStates = model.getInitialStates();
    labelingToUmb(model.getStateLabeling(), umbModel);
    if (model.hasChoiceLabeling()) {
        labelingToUmb(model.getChoiceLabeling(), umbModel);
    }
    bool normalize = model.isOfType(storm::models::ModelType::Ctmc);
    // TODO: Handle ctmc somehow
    if (!storm::NumberTraits<ValueType>::IsExact && storm::NumberTraits<TargetValueType>::IsExact) {
        STORM_LOG_WARN("Translating from non-exact to exact model representation. This may lead to rounding errors.");
        normalize = true;
    }
    transitionMatrixToUmb<ValueType, TargetValueType>(model.getTransitionMatrix(), umbModel, normalize);
    for (auto const& [name, rewardModel] : model.getRewardModels()) {
        rewardToUmb<ValueType, TargetValueType>(name, rewardModel, model.getTransitionMatrix(), umbModel);
    }
}

}  // namespace detail

template<typename ValueType>
std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb(storm::models::sparse::Model<ValueType> const& model, ExportOptions const& options) {
    auto umbModel = std::make_unique<UmbModel<StorageType::Memory>>();
    detail::setIndexInformation(model, umbModel->index, options);
    using enum ExportOptions::ValueType;
    switch (options.valueType) {
        case Default:
            detail::sparseModelToUmb<ValueType, ValueType>(model, *umbModel, options);
            break;
        case Double:
            detail::sparseModelToUmb<ValueType, double>(model, *umbModel, options);
            break;
        case Rational:
            detail::sparseModelToUmb<ValueType, storm::RationalNumber>(model, *umbModel, options);
            break;
        case DoubleInterval:
            detail::sparseModelToUmb<ValueType, storm::Interval>(model, *umbModel, options);
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unexpected value type.");
    }
    STORM_LOG_ASSERT(umbModel->validate(), "Created umb model is not valid.");
    return std::make_unique<UmbModelBase>(std::move(umbModel));
}

template std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb<double>(storm::models::sparse::Model<double> const& model, ExportOptions const& options);
template std::unique_ptr<storm::umb::UmbModelBase> sparseModelToUmb<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                                           ExportOptions const& options);
}  // namespace storm::umb
#include "storm/storage/umb/export/SparseModelToUmb.h"

#include "storm/storage/SparseMatrix.h"
#include "storm/storage/umb/model/StringEncoding.h"
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
#include "storm/storage/sparse/ChoiceOrigins.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

namespace storm::umb {

namespace detail {

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

void stateLabelingToUmb(storm::models::sparse::StateLabeling const& labeling, storm::umb::UmbModel<StorageType::Memory>& umb) {
    for (auto const& labelName : labeling.getLabels()) {
        if (labelName == "init") {
            continue;  // skip initial state labeling. Initial states are handled separately.
        }
        auto const name = umb.index.annotations.findAtomicPropositionName(labelName);
        STORM_LOG_ASSERT(name.has_value(), "Label '" << labelName << "' not found in the model index.");
        auto& annotation = umb.aps[*name];
        annotation.forStates.emplace().values.template set<bool>(labeling.getStates(labelName));
    }
}

void choiceOriginsToUmb(storm::storage::sparse::ChoiceOrigins const& choiceOrigins, storm::umb::UmbModel<StorageType::Memory>& umb) {
    // choice to action
    auto& choiceToAction = umb.choices.choiceToAction.emplace();
    choiceToAction.reserve(choiceOrigins.getNumberOfChoices());
    for (uint64_t c = 0; c < choiceOrigins.getNumberOfChoices(); ++c) {
        auto const& choiceId = choiceOrigins.getIdentifier(c);
        choiceToAction.push_back(choiceOrigins.getIdentifier(c));
    }

    // action strings
    // We always set a csr, even in cases where it could be omitted.
    // We use the empty action string for choices with no origin, which we add to the action strings by initializing the csr with {0,0}
    auto& chars = umb.choices.actionStrings.emplace();
    auto& csr = umb.choices.actionToActionString.emplace(2, 0);
    STORM_LOG_ASSERT(choiceOrigins.getIdentifierForChoicesWithNoOrigin() == 0, "Identifier for choices with no origin expected to be 0");

    for (uint64_t id = 1; id < choiceOrigins.getNumberOfIdentifiers(); ++id) {  // intentionally start at 1, since we already added the empty action string
        auto const& actionString = choiceOrigins.getIdentifierInfo(id);
        chars.insert(chars.end(), actionString.begin(), actionString.end());
        csr.push_back(chars.size());
    }
}

void choiceLabelingToUmb(storm::models::sparse::ChoiceLabeling const& labeling, storm::umb::UmbModel<StorageType::Memory>& umb) {
    // initialize umb data
    auto& choiceToAction = umb.choices.choiceToAction.emplace(labeling.getNumberOfItems(), 0);  // by default, all choices have action id 0
    auto& chars = umb.choices.actionStrings.emplace();
    auto& csr = umb.choices.actionToActionString.emplace(1, 0);  // csr mapping must start with 0. We always set a csr, even in cases where it could be omitted.

    // Auxiliary function to get the action index for a given action name and fill
    using ActionIndexType = std::remove_cvref_t<decltype(choiceToAction.front())>;
    auto getActionIndex = [&chars, &csr](std::string_view actionName) -> ActionIndexType {
        // check if action is already known or append new one
        auto actionStrings = storm::umb::stringVectorView(chars, csr);
        ActionIndexType i = 0;
        for (auto const& actionString : actionStrings) {
            if (actionString == actionName) {
                return i;  // action already exists
            }
            ++i;
        }
        // action does not exist, add it
        chars.insert(chars.end(), actionName.begin(), actionName.end());
        csr.push_back(chars.size());
        return actionStrings.size();
    };

    // Find out which choices have zero, at least one, or multiple labels. The former two cases can be handled more efficiently
    auto const labels = labeling.getLabels();
    storm::storage::BitVector choicesWithAtLeastOneLabel, choicesWithMultipleLabels;
    for (auto const& labelName : labels) {
        auto const& currentChoices = labeling.getChoices(labelName);
        if (choicesWithAtLeastOneLabel.size() == 0) {
            // first processed label
            choicesWithAtLeastOneLabel = currentChoices;
        } else if (choicesWithMultipleLabels.size() == 0) {
            // second processed label
            choicesWithMultipleLabels = choicesWithAtLeastOneLabel & currentChoices;
            choicesWithAtLeastOneLabel |= currentChoices;
        } else {
            // third or later processed label
            choicesWithMultipleLabels |= choicesWithAtLeastOneLabel & currentChoices;
            choicesWithAtLeastOneLabel |= currentChoices;
        }
    }

    // Handle choices without any labels.
    if (!choicesWithAtLeastOneLabel.full()) {
        // For consistency, unlabelled choices shall always have action index 0. So we add the empty action string.
        getActionIndex("");
        STORM_LOG_ASSERT(getActionIndex("") == 0, "Action index for empty action string must be 0.");
        // nothing else to do for unlabeled choices: we already initialized the choiceToAction mapping with 0s
    }

    // Handle choices with exactly one label.
    auto setChoices = [&choiceToAction, &getActionIndex](storm::storage::BitVector const& choices, std::string_view actionName) {
        auto choiceIt = choices.begin();
        auto const choiceItEnd = choices.end();
        if (choiceIt != choiceItEnd) {
            // there is at least one choice with this label
            auto const actionIndex = getActionIndex(actionName);
            for (; choiceIt != choiceItEnd; ++choiceIt) {
                choiceToAction[*choiceIt] = actionIndex;  // set action index for this choice
            }
        }
    };
    if (choicesWithMultipleLabels.empty()) {
        for (auto const& labelName : labels) {
            setChoices(labeling.getChoices(labelName), labelName);
        }
    } else {
        choicesWithMultipleLabels.complement();  // now contains the choices with at most one label
        for (auto const& labelName : labels) {
            setChoices(labeling.getChoices(labelName) & choicesWithMultipleLabels, labelName);
        }
        choicesWithMultipleLabels.complement();  // revert above complement operation
    }

    // Handle choices with multiple labels.
    for (auto const& choice : choicesWithMultipleLabels) {
        std::string action;
        for (auto const& label : labeling.getLabelsOfChoice(choice)) {
            if (!action.empty()) {
                action += ",";  // separate multiple labels with a comma
            }
            action += label;
        }
        choiceToAction[choice] = getActionIndex(action);
    }
}

template<typename TargetValueType>
void setGenericVector(storm::umb::GenericVector<storm::umb::StorageType::Memory>& target, std::ranges::input_range auto&& values) {
    using ValueType = std::ranges::range_value_t<decltype(values)>;
    if constexpr (std::is_same_v<ValueType, TargetValueType>) {
        target.template set<TargetValueType>(std::forward<decltype(values)>(values));
    } else {
        target.template set<TargetValueType>(storm::utility::vector::convertNumericVector<TargetValueType>(std::forward<decltype(values)>(values)));
    }
}

template<typename ValueType, typename TargetValueType>
void rewardToUmb(std::string const& identifier, storm::models::sparse::StandardRewardModel<ValueType> const& rewardModel,
                 storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::umb::UmbModel<StorageType::Memory>& umb) {
    auto const rewardName = umb.index.annotations.findRewardName(identifier);
    STORM_LOG_ASSERT(rewardName.has_value(), "Reward '" << identifier << "' not found in the model index.");
    // auto const& rewardIndex = umb.index.annotations.rewards.at(rewardName);
    STORM_LOG_ASSERT(!umb.rewards.contains(*rewardName), "Reward '" << identifier << "' already exists in the umb model.");
    auto& rewardAnnotation = umb.rewards[*rewardName];
    if (rewardModel.hasStateRewards()) {
        setGenericVector<TargetValueType>(rewardAnnotation.forStates.emplace().values, rewardModel.getStateRewardVector());
    }
    if (rewardModel.hasStateActionRewards()) {
        setGenericVector<TargetValueType>(rewardAnnotation.forChoices.emplace().values, rewardModel.getStateActionRewardVector());
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
    ts.numActions =
        (model.hasChoiceLabeling() || model.hasChoiceOrigins()) ? ts.InvalidNumber : 1;  // action count is only known after processing choice labeling/origins.
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
    using Rew = storm::umb::ModelIndex::Annotations::Annotation::Type;
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
    if (ts.time != storm::umb::ModelIndex::TransitionSystem::Time::Discrete) {
        ts.exitRateType = ts.branchProbabilityType;
    }

    // annotations:
    // rewards:
    if (model.hasRewardModel()) {
        auto& rewards = index.annotations.rewards.emplace();
        for (auto const& [rewardModelName, rewardModel] : model.getRewardModels()) {
            auto identifier = umb::ModelIndex::Annotations::Annotation::getValidIdentifierFromAlias(rewardModelName);
            STORM_LOG_THROW(!rewards.contains(identifier), storm::exceptions::WrongFormatException, "Reward id '" << identifier << "' already exists.");
            auto& rewardIndex = rewards[identifier];
            if (!rewardModelName.empty()) {
                rewardIndex.alias = rewardModelName;  // Don't introduce an alias for unnamed rewards. They don't have a nice name.
            }
            rewardIndex.type = rewardType;
            rewardIndex.appliesTo.emplace();
            using enum storm::umb::ModelIndex::Annotations::Annotation::AppliesTo;
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
    }

    // aps:
    if (model.getStateLabeling().getNumberOfLabels() > 0 || model.hasChoiceLabeling()) {
        auto& aps = index.annotations.aps.emplace();
        for (auto const& label : model.getStateLabeling().getLabels()) {
            if (label == "init") {
                continue;
            }
            auto identifier = umb::ModelIndex::Annotations::Annotation::getValidIdentifierFromAlias(label);
            STORM_LOG_THROW(!aps.contains(identifier), storm::exceptions::WrongFormatException, "AP with identifier '" << identifier << "' already exists.");
            auto& apIndex = aps[identifier];
            apIndex.alias = label;
            apIndex.type = storm::umb::ModelIndex::Annotations::Annotation::Type::Bool;
            apIndex.appliesTo.emplace();
            apIndex.appliesTo->push_back(storm::umb::ModelIndex::Annotations::Annotation::AppliesTo::States);
        }
    }
}

template<typename ValueType, typename TargetValueType>
void sparseModelToUmb(storm::models::sparse::Model<ValueType> const& model, UmbModel<StorageType::Memory>& umbModel, ExportOptions const& options) {
    setIndexInformation(model, umbModel.index, options);
    umbModel.states.initialStates = model.getInitialStates();
    stateLabelingToUmb(model.getStateLabeling(), umbModel);
    if (options.allowChoiceOriginsAsActions && model.hasChoiceOrigins()) {
        STORM_LOG_WARN_COND(!options.allowChoiceLabelingAsActions || !model.hasChoiceLabeling(),
                            "Choice origins and choice labeling are both present but only choice origins will be used as actions for UMB export.");
        choiceOriginsToUmb(*model.getChoiceOrigins(), umbModel);
        umbModel.index.transitionSystem.numActions = model.getChoiceOrigins()->getNumberOfIdentifiers();
    } else if (options.allowChoiceLabelingAsActions && model.hasChoiceLabeling()) {
        choiceLabelingToUmb(model.getChoiceLabeling(), umbModel);
        umbModel.index.transitionSystem.numActions = umbModel.choices.actionToActionString->size() - 1;
    }
    if (model.hasStateValuations()) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "State valuations are not yet supported for UMB export.");
    }
    using enum storm::models::ModelType;
    bool normalize = model.isOfType(Ctmc);
    if (!storm::NumberTraits<ValueType>::IsExact && storm::NumberTraits<TargetValueType>::IsExact) {
        STORM_LOG_WARN("Translating from non-exact to exact model representation. This may lead to rounding errors.");
        normalize = true;
    }
    transitionMatrixToUmb<ValueType, TargetValueType>(model.getTransitionMatrix(), umbModel, normalize);
    for (auto const& [name, rewardModel] : model.getRewardModels()) {
        rewardToUmb<ValueType, TargetValueType>(name, rewardModel, model.getTransitionMatrix(), umbModel);
    }
    if (model.isOfType(Ctmc)) {
        auto const& ctmc = *model.template as<storm::models::sparse::Ctmc<ValueType>>();
        setGenericVector<TargetValueType>(umbModel.states.exitRates, ctmc.getExitRateVector());
    } else if (model.isOfType(MarkovAutomaton)) {
        auto const& ma = *model.template as<storm::models::sparse::MarkovAutomaton<ValueType>>();
        umbModel.states.markovianStates = ma.getMarkovianStates();
        setGenericVector<TargetValueType>(umbModel.states.exitRates, ma.getExitRates());
    } else if (model.isOfType(Smg)) {
        auto const& smg = *model.template as<storm::models::sparse::Smg<ValueType>>();
        auto const& playerIds = smg.getStatePlayerIndications();
        umbModel.states.stateToPlayer.emplace(playerIds.begin(), playerIds.end());
    } else {
        STORM_LOG_THROW(model.isOfType(Dtmc) || model.isOfType(Mdp), storm::exceptions::NotSupportedException,
                        "Unexpected model type for UMB export: " << model.getType());
    }
}

}  // namespace detail

template<typename ValueType>
storm::umb::UmbModelBase sparseModelToUmb(storm::models::sparse::Model<ValueType> const& model, ExportOptions const& options) {
    auto umbModel = std::make_unique<UmbModel<StorageType::Memory>>();
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
    return storm::umb::UmbModelBase(std::move(umbModel));
}

template storm::umb::UmbModelBase sparseModelToUmb<double>(storm::models::sparse::Model<double> const& model, ExportOptions const& options);
template storm::umb::UmbModelBase sparseModelToUmb<storm::RationalNumber>(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                          ExportOptions const& options);
}  // namespace storm::umb
#include "storm/storage/umb/import/SparseModelFromUmb.h"

#include <ranges>
#include <utility>

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
#include "storm/storage/umb/model/StringEncoding.h"
#include "storm/storage/umb/model/ValueEncoding.h"
#include "storm/utility/builder.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/WrongFormatException.h"

namespace storm::umb {

namespace detail {

auto csrRange(auto&& csr, uint64_t i) {
    if (csr) {
        STORM_LOG_ASSERT(i + 1 < csr->size(), "CSR index out of bounds: " << (i + 1) << " >= " << csr->size());
        return std::ranges::iota_view(csr.value()[i], csr.value()[i + 1]);
    } else {
        // assume 1:1 mapping
        return std::ranges::iota_view(i, i + 1);
    }
}

template<typename ValueType, StorageType Storage>
storm::storage::SparseMatrix<ValueType> createMatrix(storm::umb::UmbModel<Storage> const& umbModel, std::ranges::input_range auto&& branchValues) {
    auto const& tsIndex = umbModel.index.transitionSystem;
    bool const hasRowGroups = tsIndex.numPlayers >= 1;
    storm::storage::SparseMatrixBuilder<ValueType> builder(tsIndex.numChoices, tsIndex.numStates, tsIndex.numBranches, true, hasRowGroups,
                                                           hasRowGroups ? tsIndex.numStates : 0u);
    for (uint64_t stateIndex{0}; stateIndex < tsIndex.numStates; ++stateIndex) {
        auto choices = csrRange(umbModel.states.stateToChoice, stateIndex);
        if (hasRowGroups) {
            builder.newRowGroup(*choices.begin());
        }
        for (auto const choiceIndex : choices) {
            STORM_LOG_ASSERT(choiceIndex < tsIndex.numChoices, "Choice index out of bounds.");
            for (auto const branchIndex : csrRange(umbModel.choices.choiceToBranch, choiceIndex)) {
                auto const& branchTarget = umbModel.branches.branchToTarget.value()[branchIndex];
                STORM_LOG_ASSERT(branchTarget < tsIndex.numStates, "Branch target index out of bounds: " << branchTarget << " >= " << tsIndex.numStates);
                builder.addNextValue(choiceIndex, branchTarget, branchValues[branchIndex]);
            }
        }
    }
    return builder.build();
}

template<typename ValueType, StorageType Storage>
storm::storage::SparseMatrix<ValueType> createMatrix(storm::umb::UmbModel<Storage> const& umbModel, auto sourceType,
                                                     storm::umb::GenericVector<Storage> const& branchValues,
                                                     typename storm::umb::UmbModel<Storage>::CSR const& csr) {
    return ValueEncoding::applyDecodedVector<ValueType, Storage>([&umbModel](auto&& input) { return createMatrix<ValueType, Storage>(umbModel, input); },
                                                                 branchValues, sourceType, csr);
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
storm::storage::BitVector createBitVector(std::optional<storm::umb::VectorType<bool, Storage>> const& umbBitVector, uint64_t size) {
    STORM_LOG_THROW(umbBitVector.has_value(), storm::exceptions::WrongFormatException, "BitVector is not given but expected.");
    return createBitVector<Storage>(*umbBitVector, size);
}

template<StorageType Storage>
storm::models::sparse::StateLabeling constructStateLabeling(storm::umb::UmbModel<Storage> const& umbModel) {
    auto const& numStates = umbModel.index.transitionSystem.numStates;
    storm::models::sparse::StateLabeling stateLabelling(numStates);
    if (umbModel.states.initialStates) {
        stateLabelling.addLabel("init", createBitVector<Storage>(umbModel.states.initialStates, numStates));
    } else {
        STORM_LOG_WARN("No initial states given in UMB model.");
        stateLabelling.addLabel("init", storm::storage::BitVector(numStates, false));  // default to all states not being initial
    }
    for (auto const& [apName, ap] : umbModel.aps) {
        STORM_LOG_THROW(umbModel.index.annotations.aps->contains(apName), storm::exceptions::WrongFormatException,
                        "Atomic proposition '" << apName << "' not found in index.");
        auto const& apIndex = umbModel.index.annotations.aps->at(apName);
        auto labelName = apIndex.alias.value_or(apName);  // prefer alias as label name if it exists
        STORM_LOG_THROW(ap.forStates.has_value(), storm::exceptions::WrongFormatException, "Atomic proposition '" << apName << "' does not apply to states.");
        STORM_LOG_THROW(!stateLabelling.containsLabel(labelName), storm::exceptions::WrongFormatException,
                        "Label '" << labelName << "' already exists in state labeling.");
        stateLabelling.addLabel(labelName, createBitVector<Storage>(ap.forStates->values.template get<bool>(), numStates));
        STORM_LOG_WARN_COND(!ap.forChoices.has_value(), "Atomic propositions for choices are not supported.");
        STORM_LOG_WARN_COND(!ap.forBranches.has_value(), "Atomic propositions for branches are not supported.");
    }
    return stateLabelling;
}

template<StorageType Storage>
storm::models::sparse::ChoiceLabeling constructChoiceLabeling(storm::umb::UmbModel<Storage> const& umbModel) {
    auto const& numChoices = umbModel.index.transitionSystem.numChoices;
    storm::models::sparse::ChoiceLabeling choiceLabeling(numChoices);
    auto actionStrings = storm::umb::stringVectorViewOptionalCsr(umbModel.choices.actionStrings.value(), umbModel.choices.actionToActionString);

    uint64_t const emptyActionIndex = std::ranges::find(actionStrings, std::string_view{}) - actionStrings.begin();

    // for each choice, find the corresponding action index and set the bit accordingly
    auto const& choiceToAction = umbModel.choices.choiceToAction.value();
    std::vector<storm::storage::BitVector> actionToLabels(actionStrings.size(), storm::storage::BitVector(numChoices, false));
    for (uint64_t choiceIndex = 0; choiceIndex < numChoices; ++choiceIndex) {
        auto const actionIndex = choiceToAction[choiceIndex];
        STORM_LOG_ASSERT(actionIndex < actionStrings.size(), "Choice to action mapping out of bounds.");
        if (actionIndex == emptyActionIndex) {
            continue;  // skip choices with empty action. They will not be labeled.
        }
        actionToLabels[actionIndex].set(choiceIndex);
    }

    // add the action labels to the labeling
    for (uint64_t actionIndex = 0; actionIndex < actionStrings.size(); ++actionIndex) {
        if (actionIndex == emptyActionIndex) {
            continue;
        }
        choiceLabeling.addLabel(std::string(actionStrings[actionIndex]), std::move(actionToLabels[actionIndex]));
    }
    return choiceLabeling;
}

template<typename ValueType, StorageType Storage>
auto constructRewardModels(storm::umb::UmbModel<Storage> const& umbModel) {
    using RewardModel = storm::models::sparse::StandardRewardModel<ValueType>;
    std::unordered_map<std::string, RewardModel> rewardModels;
    for (auto const& [rewName, rew] : umbModel.rewards) {
        STORM_LOG_THROW(umbModel.rewards.contains(rewName), storm::exceptions::WrongFormatException, "Reward '" << rewName << "' not found in index.");
        auto const& rewIndex = umbModel.index.annotations.rewards->at(rewName);
        auto usedRewName = rewIndex.alias.value_or(rewName);  // prefer alias as reward name if it exists
        STORM_LOG_THROW(!rewardModels.contains(usedRewName), storm::exceptions::WrongFormatException,
                        "Reward '" << usedRewName << "' already exists in reward models.");
        STORM_LOG_THROW(rewIndex.type.has_value(), storm::exceptions::WrongFormatException,
                        "Reward type is not set in the index for reward '" << rewName << "'.");
        std::optional<std::vector<ValueType>> stateRewards, stateActionRewards;
        std::optional<storm::storage::SparseMatrix<ValueType>> transitionRewards;
        if (rew.forStates) {
            stateRewards = ValueEncoding::createDecodedVector<ValueType, Storage>(rew.forStates->values, rewIndex.type.value(), rew.forStates->toValue);
        }
        if (rew.forChoices) {
            stateActionRewards = ValueEncoding::createDecodedVector<ValueType, Storage>(rew.forChoices->values, rewIndex.type.value(), rew.forChoices->toValue);
        }
        if (rew.forBranches) {
            transitionRewards = createMatrix<ValueType, Storage>(umbModel, rewIndex.type.value(), rew.forBranches->values, rew.forBranches->toValue);
        }
        rewardModels.emplace(std::move(usedRewName), RewardModel(std::move(stateRewards), std::move(stateActionRewards), std::move(transitionRewards)));
    }
    return rewardModels;
}

template<typename ValueType, StorageType Storage>
std::shared_ptr<storm::models::sparse::Model<ValueType>> constructSparseModel(storm::umb::UmbModel<Storage> const& umbModel, ImportOptions const& options) {
    STORM_LOG_THROW(umbModel.validate(), storm::exceptions::WrongFormatException, "UMB model is not valid.");

    // transitions, labelings, rewards
    auto stateLabelling = constructStateLabeling(umbModel);
    auto transitionMatrix = createMatrix<ValueType>(umbModel, umbModel.index.transitionSystem.branchProbabilityType, umbModel.branches.branchProbabilities,
                                                    umbModel.branches.branchToProbability);
    storm::storage::sparse::ModelComponents<ValueType> components(std::move(transitionMatrix), std::move(stateLabelling),
                                                                  constructRewardModels<ValueType>(umbModel));
    if (options.buildChoiceLabeling && umbModel.choices.choiceToAction) {
        components.choiceLabeling = constructChoiceLabeling(umbModel);
    }
    if (options.buildStateValuations) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "State valuations are not yet supported.");
    }

    // model type-specific components
    using enum storm::models::ModelType;
    auto const modelType = deriveModelType(umbModel.index);
    if (modelType == Ctmc || modelType == MarkovAutomaton) {
        STORM_LOG_THROW(umbModel.states.exitRates.hasValue(), storm::exceptions::WrongFormatException,
                        "Exit rates are required for CTMC and Markov automaton models but not present in the UMB model.");
        components.exitRates = ValueEncoding::createDecodedVector<ValueType, Storage>(
            umbModel.states.exitRates, umbModel.index.transitionSystem.exitRateType.value(), umbModel.states.stateToExitRate);
        if (modelType == MarkovAutomaton) {
            STORM_LOG_THROW(umbModel.states.markovianStates.has_value(), storm::exceptions::WrongFormatException,
                            "Markovian states are required for Markov automaton models but not present in the UMB model.");
            components.markovianStates = createBitVector<Storage>(umbModel.states.markovianStates, umbModel.index.transitionSystem.numStates);
        }
    } else if (modelType == Smg) {
        STORM_LOG_THROW(umbModel.states.stateToPlayer.has_value(), storm::exceptions::WrongFormatException,
                        "Player information is required for SMG models but not present in the UMB model.");
        auto const stateToPlayer = umbModel.states.stateToPlayer.value();
        components.statePlayerIndications.emplace(stateToPlayer.begin(), stateToPlayer.end());
    } else {
        STORM_LOG_THROW(modelType == Dtmc || modelType == Mdp, storm::exceptions::NotSupportedException, "Unexpected model type for UMB import: " << modelType);
    }
    return storm::utility::builder::buildModelFromComponents(deriveModelType(umbModel.index), std::move(components));
}

}  // namespace detail

storm::models::ModelType deriveModelType(storm::umb::ModelIndex const& index) {
    using enum storm::models::ModelType;

    auto const& ts = index.transitionSystem;

    STORM_LOG_THROW(ts.branchProbabilityType != storm::umb::ModelIndex::TransitionSystem::BranchProbabilityType::None, storm::exceptions::NotSupportedException,
                    "Models without branch values are not supported.");
    switch (ts.time) {
        using enum storm::umb::ModelIndex::TransitionSystem::Time;
        case Discrete:
            switch (ts.numPlayers) {
                case 0:
                    return Dtmc;
                case 1:
                    return Mdp;
                default:
                    return Smg;
            }
        case Stochastic:
            STORM_LOG_THROW(ts.numPlayers == 0, storm::exceptions::NotSupportedException, "Stochastic time models with multiple players are not supported.");
            return Ctmc;
        case UrgentStochastic:
            STORM_LOG_THROW(ts.numPlayers == 1, storm::exceptions::NotSupportedException,
                            "Urgent stochastic time models with multiple or no players are not supported.");
            return MarkovAutomaton;
    }
    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected transition system time type" << ts.time << ".");
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModelFromUmb(storm::umb::UmbModelBase const& umbModel, ImportOptions const& options) {
    if (umbModel.isStorageType(StorageType::Disk)) {
        return detail::constructSparseModel<ValueType, StorageType::Disk>(umbModel.as<StorageType::Disk>(), options);
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        return detail::constructSparseModel<ValueType, StorageType::Memory>(umbModel.as<StorageType::Memory>(), options);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Storage type not expected.");
    }
}
template std::shared_ptr<storm::models::sparse::Model<double>> sparseModelFromUmb<double>(storm::umb::UmbModelBase const& umbModel,
                                                                                          ImportOptions const& options);
template std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> sparseModelFromUmb<storm::RationalNumber>(
    storm::umb::UmbModelBase const& umbModel, ImportOptions const& options);
}  // namespace storm::umb
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
#include "storm/utility/builder.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/WrongFormatException.h"

namespace storm::umb {

namespace detail {

template<typename TargetType, typename SourceType>
auto conversionView(std::ranges::input_range auto&& input) {
    if constexpr (std::same_as<TargetType, SourceType>) {
        return input;
    } else {
        return input | std::ranges::views::transform(
                           [](SourceType const& value) -> TargetType { return storm::utility::convertNumber<TargetType, SourceType>(value); });
    }
}

template<StorageType Storage>
bool hasType(storm::umb::GenericVector<Storage> const& input, auto type) {
    switch (type) {
        case decltype(type)::Double:
            return input.template isType<double>();
        case decltype(type)::Rational:
            return input.template isType<storm::RationalNumber>();
        case decltype(type)::DoubleInterval:
            return input.template isType<storm::Interval>();
        default:
            return false;
    }
}

template<auto Source, bool IsInSourceTypeRepresentation, typename TargetType, StorageType Storage>
auto valueVectorView(storm::umb::GenericVector<Storage> const& input, typename storm::umb::UmbModel<Storage>::CSR const& csr = {}) {
    STORM_LOG_ASSERT(IsInSourceTypeRepresentation == hasType(input, Source), "Inconsistent arguments.");
    if constexpr (Source == decltype(Source)::Double) {
        static_assert(IsInSourceTypeRepresentation, "Double is the only valid source type in this case.");
        STORM_LOG_ASSERT(input.template isType<double>(), "Unexpected type for values. Expected double.");
        STORM_LOG_WARN_COND(!storm::NumberTraits<TargetType>::IsExact,
                            "Some values are given in type double but will be converted to an exact (arbitrary precision) type. Rounding errors may occur.");
        return conversionView<TargetType, double>(input.template get<double>());

    } else if constexpr (Source == decltype(Source)::Rational) {
        STORM_LOG_WARN_COND(storm::NumberTraits<TargetType>::IsExact,
                            "Some values are given in an exact type but converted to an inexact type. Rounding errors may occur.");
        if constexpr (IsInSourceTypeRepresentation) {
            return conversionView<TargetType, storm::RationalNumber>(input.template get<storm::RationalNumber>());
        } else {
            STORM_LOG_ASSERT(input.template isType<uint64_t>(), "Unexpected type for rational representation. Expected uint64.");
            // TODO: complete
            assert(false);
            return std::vector<TargetType>{};
        }
    } else if constexpr (Source == decltype(Source)::DoubleInterval) {
        if constexpr (!std::is_same_v<TargetType, storm::Interval>) {
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException,
                            "Some values are given as double intervals but a model with a non-interval type is requested.");
            return std::vector<TargetType>{};
        } else if constexpr (IsInSourceTypeRepresentation) {
            return input.template get<storm::Interval>();
        } else {
            STORM_LOG_ASSERT(input.template isType<double>(), "Unexpected type for double interval representation. Expected double.");
            // TODO: complete
            assert(false);
            return std::vector<TargetType>{};
        }
    } else {
        static_assert(false, "Unexpected type for values.");
    }
}

template<typename ValueType, StorageType Storage>
auto createHelper(auto sourceType, storm::umb::GenericVector<Storage> const& input, typename storm::umb::UmbModel<Storage>::CSR const& csr, auto&& create) {
    using E = decltype(sourceType)::E;

    // find out how to interpret the input values
    switch (sourceType) {
        case E::Double:
            return create(valueVectorView<E::Double, true, ValueType, Storage>(input));
        case E::Rational:
            if (hasType(input, E::Rational)) {
                return create(valueVectorView<E::Rational, true, ValueType, Storage>(input));
            } else {
                // Only this case might require the csr. It is useless in all other cases.
                return create(valueVectorView<E::Rational, false, ValueType, Storage>(input, csr));
            }
        case E::DoubleInterval:
            if (hasType(input, E::DoubleInterval)) {
                return create(valueVectorView<E::DoubleInterval, true, ValueType, Storage>(input));
            } else {
                return create(valueVectorView<E::DoubleInterval, false, ValueType, Storage>(input));
            }
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Values have unsupported type " << sourceType << ".");
    }
}

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
    return createHelper<ValueType, Storage>(sourceType, branchValues, csr,
                                            [&umbModel](auto&& input) { return createMatrix<ValueType, Storage>(umbModel, input); });
}

template<typename ValueType, StorageType Storage>
std::vector<ValueType> createVector(auto sourceType, storm::umb::GenericVector<Storage> const& input, typename storm::umb::UmbModel<Storage>::CSR const& csr) {
    return createHelper<ValueType, Storage>(sourceType, input, csr, [](auto&& input) { return std::vector<ValueType>(input.begin(), input.end()); });
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
std::pair<storm::models::sparse::StateLabeling, std::optional<storm::models::sparse::ChoiceLabeling>> constructStateChoiceLabelling(
    storm::umb::UmbModel<Storage> const& umbModel) {
    auto const& tsIndex = umbModel.index.transitionSystem;
    storm::models::sparse::StateLabeling stateLabelling(tsIndex.numStates);
    std::optional<storm::models::sparse::StateLabeling> choiceLabelling;
    if (umbModel.states.initialStates) {
        stateLabelling.addLabel("init", createBitVector<Storage>(umbModel.states.initialStates, tsIndex.numStates));
    } else {
        STORM_LOG_WARN("No initial states given in UMB model.");
    }
    for (auto const& [apName, ap] : umbModel.atomicPropositions) {
        STORM_LOG_THROW(umbModel.index.annotations.atomicPropositions.contains(apName), storm::exceptions::WrongFormatException,
                        "Atomic proposition '" << apName << "' not found in index.");
        auto const& apIndex = umbModel.index.annotations.atomicPropositions.at(apName);
        auto labelName = apIndex.alias.value_or(apName);  // prefer alias as label name if it exists
        if (ap.forStates) {
            STORM_LOG_THROW(!stateLabelling.containsLabel(labelName), storm::exceptions::WrongFormatException,
                            "Label '" << labelName << "' already exists in state labeling.");
            stateLabelling.addLabel(labelName, createBitVector<Storage>(ap.forStates->values, tsIndex.numStates));
        }
        if (ap.forChoices) {
            if (!choiceLabelling.has_value()) {
                choiceLabelling.emplace(tsIndex.numChoices);
            }
            STORM_LOG_THROW(!choiceLabelling->containsLabel(labelName), storm::exceptions::WrongFormatException,
                            "Label '" << labelName << "' already exists in choice labeling.");
            choiceLabelling->addLabel(labelName, createBitVector<Storage>(ap.forChoices->values, tsIndex.numChoices));
        }
        STORM_LOG_WARN_COND(!ap.forBranches.has_value(), "Atomic propositions for branches are not supported.");
    }
    return {std::move(stateLabelling), std::move(choiceLabelling)};
}

template<typename ValueType, StorageType Storage>
auto constructRewardModels(storm::umb::UmbModel<Storage> const& umbModel) {
    using RewardModel = storm::models::sparse::StandardRewardModel<ValueType>;
    auto const& ts = umbModel.index.transitionSystem;
    std::unordered_map<std::string, RewardModel> rewardModels;
    for (auto const& [rewName, rew] : umbModel.rewards) {
        STORM_LOG_THROW(umbModel.rewards.contains(rewName), storm::exceptions::WrongFormatException, "Reward '" << rewName << "' not found in index.");
        auto const& rewIndex = umbModel.index.annotations.rewards.at(rewName);
        auto usedRewName = rewIndex.alias.value_or(rewName);  // prefer alias as reward name if it exists
        STORM_LOG_THROW(!rewardModels.contains(usedRewName), storm::exceptions::WrongFormatException,
                        "Reward '" << usedRewName << "' already exists in reward models.");
        std::optional<std::vector<ValueType>> stateRewards, stateActionRewards;
        std::optional<storm::storage::SparseMatrix<ValueType>> transitionRewards;
        if (rew.forStates) {
            stateRewards = createVector<ValueType, Storage>(rewIndex.type, rew.forStates->values, rew.forStates->toValue);
        }
        if (rew.forChoices) {
            stateActionRewards = createVector<ValueType, Storage>(rewIndex.type, rew.forChoices->values, rew.forChoices->toValue);
        }
        if (rew.forBranches) {
            transitionRewards = createMatrix<ValueType, Storage>(umbModel, rewIndex.type, rew.forBranches->values, rew.forBranches->toValue);
        }
        rewardModels.emplace(std::move(usedRewName), RewardModel(std::move(stateRewards), std::move(stateActionRewards), std::move(transitionRewards)));
    }
    return rewardModels;
}

template<typename ValueType, StorageType Storage>
std::shared_ptr<storm::models::sparse::Model<ValueType>> constructSparseModel(storm::umb::UmbModel<Storage> const& umbModel) {
    STORM_LOG_THROW(umbModel.validate(), storm::exceptions::WrongFormatException, "UMB model is not valid.");
    auto [stateLabelling, choiceLabelling] = constructStateChoiceLabelling<Storage>(umbModel);
    auto transitionMatrix = createMatrix<ValueType>(umbModel, umbModel.index.transitionSystem.branchProbabilityType, umbModel.branches.branchProbabilities,
                                                    umbModel.branches.branchToProbability);
    storm::storage::sparse::ModelComponents<ValueType> components(std::move(transitionMatrix), std::move(stateLabelling),
                                                                  constructRewardModels<ValueType>(umbModel));
    components.choiceLabeling = std::move(choiceLabelling);
    return storm::utility::builder::buildModelFromComponents(deriveModelType(umbModel.index), std::move(components));
}

}  // namespace detail

storm::models::ModelType deriveModelType(storm::umb::ModelIndex const& index) {
    using ModelType = storm::models::ModelType;

    auto const& ts = index.transitionSystem;

    STORM_LOG_THROW(ts.branchProbabilityType != storm::umb::ModelIndex::TransitionSystem::BranchProbabilityType::None, storm::exceptions::NotSupportedException,
                    "Models without branch values are not supported.");
    switch (ts.time) {
        using enum storm::umb::ModelIndex::TransitionSystem::Time;
        case Discrete:
            switch (ts.numPlayers) {
                case 0:
                    return ModelType::Dtmc;
                case 1:
                    return ModelType::Mdp;
                default:
                    return ModelType::Smg;
            }
        case Stochastic:
            STORM_LOG_THROW(ts.numPlayers == 0, storm::exceptions::NotSupportedException, "Stochastic time models with multiple players are not supported.");
            return ModelType::Ctmc;
        case UrgentStochastic:
            STORM_LOG_THROW(ts.numPlayers == 1, storm::exceptions::NotSupportedException,
                            "Urgent stochastic time models with multiple or no players are not supported.");
            return ModelType::MarkovAutomaton;
    }
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModelFromUmb(storm::umb::UmbModelBase const& umbModel) {
    if (umbModel.isStorageType(StorageType::Disk)) {
        return detail::constructSparseModel<ValueType, StorageType::Disk>(umbModel.as<StorageType::Disk>());
    } else if (umbModel.isStorageType(StorageType::Memory)) {
        return detail::constructSparseModel<ValueType, StorageType::Memory>(umbModel.as<StorageType::Memory>());
    } else {
        STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Storage type not expected.");
    }
}
template std::shared_ptr<storm::models::sparse::Model<double>> sparseModelFromUmb<double>(storm::umb::UmbModelBase const& umbModel);
template std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> sparseModelFromUmb<storm::RationalNumber>(
    storm::umb::UmbModelBase const& umbModel);
}  // namespace storm::umb
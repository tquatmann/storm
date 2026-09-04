#pragma once

#include <optional>

#include "storm/exceptions/NotSupportedException.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/storage/dd/DdType.h"
#include "storm/storage/dd/bisimulation/BisimulationDecomposition.h"
#include "storm/storage/dd/bisimulation/BisimulationOptions.h"
#include "storm/storage/stateminimization/bisimulation/DeterministicModelBisimulationDecomposition.h"
#include "storm/storage/stateminimization/bisimulation/IntervalModelBisimulationDecomposition.h"
#include "storm/storage/stateminimization/bisimulation/NondeterministicModelBisimulationDecomposition.h"
#include "storm/transformer/bisimulation/Bisimulation.h"
#include "storm/utility/NumberTraits.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm {
namespace api {

template<typename ModelType>
std::shared_ptr<ModelType> performDeterministicSparseBisimulationMinimization(std::shared_ptr<ModelType> model,
                                                                              std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                              storm::storage::BisimulationType type, bool graphPreserving = true,
                                                                              std::optional<double> const& tolerance = std::nullopt) {
    // Falls back to the general precision setting when the caller does not deliberately choose a tolerance;
    // may be reworked to require an explicit choice throughout the API in the future.
    storm::IntervalBaseType<typename ModelType::ValueType> const resolvedTolerance =
        storm::NumberTraits<storm::IntervalBaseType<typename ModelType::ValueType>>::IsExact
            ? storm::utility::zero<storm::IntervalBaseType<typename ModelType::ValueType>>()
            : storm::utility::convertNumber<storm::IntervalBaseType<typename ModelType::ValueType>>(
                  tolerance.value_or(storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision()));

    // TODO: Make a clean distinction between interval and standard models here.
    // TODO: Instantiate the corresponding implementation for bisimulation.
    if constexpr (storm::IsIntervalType<typename ModelType::ValueType>) {
        using OptionsType = typename storm::storage::IntervalModelBisimulationDecomposition<ModelType>::BisimulationOptions;
        OptionsType options =
            (!formulas.empty() && graphPreserving) ? OptionsType(*model, formulas, resolvedTolerance) : OptionsType::preservingAllLabels(resolvedTolerance);
        // If we cannot use formula-based decomposition because of
        // non-graph-preserving regions but there are reward models, we need to
        // preserve those
        if (!graphPreserving &&
            std::any_of(formulas.begin(), formulas.end(), [](auto const& formula) { return formula->getReferencedRewardModels().size() > 0; })) {
            options.setKeepRewards(true);
        }
        options.setType(type);

        storm::storage::IntervalModelBisimulationDecomposition<ModelType> bisimulationDecomposition(*model, options);
        bisimulationDecomposition.computeDecomposition();
        return bisimulationDecomposition.getQuotient();
    } else {
        using OptionsType = typename storm::storage::DeterministicModelBisimulationDecomposition<ModelType>::BisimulationOptions;
        OptionsType options =
            (!formulas.empty() && graphPreserving) ? OptionsType(*model, formulas, resolvedTolerance) : OptionsType::preservingAllLabels(resolvedTolerance);
        if (!graphPreserving &&
            std::any_of(formulas.begin(), formulas.end(), [](auto const& formula) { return formula->getReferencedRewardModels().size() > 0; })) {
            options.setKeepRewards(true);
        }
        options.setType(type);

        storm::storage::DeterministicModelBisimulationDecomposition<ModelType> bisimulationDecomposition(*model, options);
        bisimulationDecomposition.computeDecomposition();
        return bisimulationDecomposition.getQuotient();
    }
}

template<typename ModelType>
std::shared_ptr<ModelType> performNondeterministicSparseBisimulationMinimization(std::shared_ptr<ModelType> model,
                                                                                 std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                                 storm::storage::BisimulationType type, bool graphPreserving = true,
                                                                                 bool actionSensitive = false,
                                                                                 std::optional<double> const& tolerance = std::nullopt) {
    // Falls back to the general precision setting when the caller does not deliberately choose a tolerance;
    // may be reworked to require an explicit choice throughout the API in the future.
    storm::IntervalBaseType<typename ModelType::ValueType> const resolvedTolerance =
        storm::NumberTraits<storm::IntervalBaseType<typename ModelType::ValueType>>::IsExact
            ? storm::utility::zero<storm::IntervalBaseType<typename ModelType::ValueType>>()
            : storm::utility::convertNumber<storm::IntervalBaseType<typename ModelType::ValueType>>(
                  tolerance.value_or(storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision()));

    if constexpr (storm::IsIntervalType<typename ModelType::ValueType>) {
        using OptionsType = typename storm::storage::IntervalModelBisimulationDecomposition<ModelType>::BisimulationOptions;
        OptionsType options =
            (!formulas.empty() && graphPreserving) ? OptionsType(*model, formulas, resolvedTolerance) : OptionsType::preservingAllLabels(resolvedTolerance);
        if (!graphPreserving &&
            std::any_of(formulas.begin(), formulas.end(), [](auto const& formula) { return formula->getReferencedRewardModels().size() > 0; })) {
            options.setKeepRewards(true);
        }
        options.setType(type);
        options.setActionSensitivity(actionSensitive);

        storm::storage::IntervalModelBisimulationDecomposition<ModelType> bisimulationDecomposition(*model, options);
        bisimulationDecomposition.computeDecomposition();
        return bisimulationDecomposition.getQuotient();
    } else {
        using OptionsType = typename storm::storage::NondeterministicModelBisimulationDecomposition<ModelType>::BisimulationOptions;
        OptionsType options =
            (!formulas.empty() && graphPreserving) ? OptionsType(*model, formulas, resolvedTolerance) : OptionsType::preservingAllLabels(resolvedTolerance);
        if (!graphPreserving &&
            std::any_of(formulas.begin(), formulas.end(), [](auto const& formula) { return formula->getReferencedRewardModels().size() > 0; })) {
            options.setKeepRewards(true);
        }
        options.setType(type);
        options.setActionSensitivity(actionSensitive);

        storm::storage::NondeterministicModelBisimulationDecomposition<ModelType> bisimulationDecomposition(*model, options);
        bisimulationDecomposition.computeDecomposition();
        return bisimulationDecomposition.getQuotient();
    }
}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> performBisimulationMinimization(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
    storm::storage::BisimulationType type = storm::storage::BisimulationType::Strong, bool graphPreserving = true, bool actionSensitive = false,
    std::optional<double> const& tolerance = std::nullopt) {
    STORM_LOG_THROW(
        model->isOfType(storm::models::ModelType::Dtmc) || model->isOfType(storm::models::ModelType::Ctmc) || model->isOfType(storm::models::ModelType::Mdp),
        storm::exceptions::NotSupportedException, "Bisimulation minimization is currently only available for DTMCs, CTMCs and MDPs.");

    if (type == storm::storage::BisimulationType::Weak) {
        std::cout << "Hack for using the new bisim impl" << std::endl;
        storm::bisimulation::Options options;
        options.tolerance = tolerance.value_or(0.0);
        std::cout << "Tolerance is " << options.tolerance << std::endl;
        return storm::bisimulation::applyBisimulationMinimization(*model, options, formulas).quotient;
    }

    // Try to get rid of non state-rewards to easy bisimulation computation.
    model->reduceToStateBasedRewards();

    if (model->isOfType(storm::models::ModelType::Dtmc)) {
        return performDeterministicSparseBisimulationMinimization<storm::models::sparse::Dtmc<ValueType>>(
            model->template as<storm::models::sparse::Dtmc<ValueType>>(), formulas, type, graphPreserving, tolerance);
    } else if (model->isOfType(storm::models::ModelType::Ctmc)) {
        return performDeterministicSparseBisimulationMinimization<storm::models::sparse::Ctmc<ValueType>>(
            model->template as<storm::models::sparse::Ctmc<ValueType>>(), formulas, type, graphPreserving, tolerance);
    } else {
        return performNondeterministicSparseBisimulationMinimization<storm::models::sparse::Mdp<ValueType>>(
            model->template as<storm::models::sparse::Mdp<ValueType>>(), formulas, type, graphPreserving, actionSensitive, tolerance);
    }
}

template<storm::dd::DdType DdType, typename ValueType, typename ExportValueType = ValueType>
typename std::enable_if<DdType == storm::dd::DdType::Sylvan || std::is_same<ValueType, double>::value,
                        std::shared_ptr<storm::models::Model<ExportValueType>>>::type
performBisimulationMinimization(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model,
                                std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                storm::storage::BisimulationType const& bisimulationType = storm::storage::BisimulationType::Strong,
                                storm::dd::bisimulation::SignatureMode const& mode = storm::dd::bisimulation::SignatureMode::Eager,
                                storm::dd::bisimulation::QuotientFormat const& quotientFormat = storm::dd::bisimulation::QuotientFormat::Dd,
                                storm::dd::bisimulation::BisimulationOptions const& bisimulationOptions = storm::dd::bisimulation::BisimulationOptions()) {
    STORM_LOG_THROW(model->isOfType(storm::models::ModelType::Dtmc) || model->isOfType(storm::models::ModelType::Ctmc) ||
                        model->isOfType(storm::models::ModelType::Mdp) || model->isOfType(storm::models::ModelType::MarkovAutomaton),
                    storm::exceptions::NotSupportedException, "Symbolic bisimulation minimization is currently only available for DTMCs, CTMCs, MDPs and MAs.");
    STORM_LOG_THROW(bisimulationType == storm::storage::BisimulationType::Strong, storm::exceptions::NotSupportedException,
                    "Currently only strong bisimulation is supported.");

    std::shared_ptr<storm::models::Model<ExportValueType>> result;
    model->getManager().execute([&]() {
        // Try to get rid of non state-rewards to easy bisimulation computation.
        model->reduceToStateBasedRewards();

        storm::dd::BisimulationDecomposition<DdType, ValueType, ExportValueType> decomposition(*model, formulas, bisimulationType, bisimulationOptions);
        decomposition.compute(mode);
        result = decomposition.getQuotient(quotientFormat);
    });
    return result;
}

template<storm::dd::DdType DdType, typename ValueType, typename ExportValueType = ValueType>
typename std::enable_if<DdType != storm::dd::DdType::Sylvan && !std::is_same<ValueType, double>::value,
                        std::shared_ptr<storm::models::Model<ExportValueType>>>::type
performBisimulationMinimization(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const&,
                                std::vector<std::shared_ptr<storm::logic::Formula const>> const&,
                                storm::storage::BisimulationType const& = storm::storage::BisimulationType::Strong,
                                storm::dd::bisimulation::SignatureMode const& = storm::dd::bisimulation::SignatureMode::Eager,
                                storm::dd::bisimulation::QuotientFormat const& = storm::dd::bisimulation::QuotientFormat::Dd,
                                storm::dd::bisimulation::BisimulationOptions const& = storm::dd::bisimulation::BisimulationOptions()) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                    "Symbolic bisimulation minimization is not supported for this combination of DD library and value type.");
    return nullptr;
}

}  // namespace api
}  // namespace storm

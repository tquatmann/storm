#include "storm/transformer/bisimulation/Bisimulation.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/models/sparse/Model.h"

#include "storm/transformer/bisimulation/Initialization.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Quotient.h"
#include "storm/transformer/bisimulation/QuotientData.h"
#include "storm/transformer/bisimulation/SignatureBasedRefinement.h"
#include "storm/transformer/bisimulation/Signatures.h"
#include "storm/transformer/bisimulation/SplitterBasedRefinement.h"
#include "storm/transformer/bisimulation/WeakBisimulationData.h"
#include "storm/utility/OptionalRef.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

template<typename ValueType, typename QuotientValueType>
ReturnType<QuotientValueType> performBisimulationMinimization(storm::models::sparse::Model<ValueType> const& model,
                                                              std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                              Options const& options) {
    // Step 0: Sanity checks and set-up
    // Interval abstraction takes a non-interval model and produces an interval quotient
    bool constexpr useIntervalAbstraction = (!storm::IsIntervalType<ValueType> && storm::IsIntervalType<QuotientValueType> &&
                                             std::is_same_v<ValueType, storm::IntervalBaseType<QuotientValueType>>);
    static_assert(std::is_same_v<ValueType, QuotientValueType> || useIntervalAbstraction,
                  "The case ValueType != QuotientValueType is only implemented for the interval-abstraction code path.");
    STORM_LOG_THROW(options.tolerance >= storm::utility::zero<storm::RationalNumber>(), storm::exceptions::InvalidArgumentException,
                    "Tolerance for bisimulation minimization must be non-negative, but was " << options.tolerance << ".");
    bool const isWeak = options.bisimulationType == BisimulationType::Weak;
    if (isWeak) {
        STORM_LOG_THROW(!model.isNondeterministicModel(), storm::exceptions::NotSupportedException,
                        "Weak bisimulation is only supported for deterministic models, but the given model is of type " << model.getType() << ".");
        STORM_LOG_WARN_COND(!options.preferSignatureRefinement,
                            "Using splitter-based refinement because weak bisimulation is not supported for signature-based refinement.");
    }
    if constexpr (storm::IsIntervalType<QuotientValueType>) {
        STORM_LOG_THROW(model.isOfType(storm::models::ModelType::Dtmc) || model.isOfType(storm::models::ModelType::Mdp),
                        storm::exceptions::NotSupportedException,
                        "Bisimulation with interval values is only supported for DTMCs and MDPs, but the given model is of type " << model.getType() << ".");
        STORM_LOG_THROW(!isWeak, storm::exceptions::NotSupportedException, "Weak bisimulation is not supported for interval values.");
    }
    // Determine whether to use signature-based refinement: Either if splitter-based is not supported or if signature-based is preferred.
    // Weak bisimulation with signature-based refinement is not supported either way.
    // Deterministic models without intervals default to splitter-based refinement, which is usually faster
    bool const useSignatureRefinement =
        !isWeak && (storm::IsIntervalType<QuotientValueType> || model.isNondeterministicModel() || options.preferSignatureRefinement);
    storm::utility::Stopwatch sw(true);
    STORM_LOG_STATISTICS("-------- Bisimulation Minimization --------");

    // Step 1: Obtain an initial partition based on what needs to be preserved (labels, rewards, ...)
    storm::bisimulation::Initialization<ValueType> initialization(model, options, formulas);
    auto const preservationInformation = initialization.getPreservationInformation();
    auto const choiceClasses = initialization.getChoiceClasses();
    auto partition = initialization.getInitialStatePartition(choiceClasses);
    STORM_LOG_STATISTICS(sw << " seconds for initial partition (" << partition.getNumberOfBlocks() << " blocks).");
    sw.restart();

    // Step 2: Apply refinement using the initial partition and choiceClasses. Initialize QuotientData.
    std::optional<storm::bisimulation::QuotientData<QuotientValueType>> quotientData;
    // commonly called right after refinement.
    auto initializeQuotientData = [&partition, &sw,
                                   &quotientData](storm::OptionalRef<storm::storage::BitVector const> preferredRepresentatives = storm::NullRef) {
        STORM_LOG_STATISTICS(sw << " seconds for refinement (" << partition.getNumberOfBlocks() << " blocks).");
        sw.restart();
        quotientData.emplace(partition, preferredRepresentatives);
    };
    if (useSignatureRefinement) {
        if (storm::utility::isZero(options.tolerance)) {
            if constexpr (!useIntervalAbstraction) {
                storm::bisimulation::Signatures<ValueType, SignatureMode::Exact, QuotientValueType> signatures(model, choiceClasses, partition);
                storm::bisimulation::performSignatureBasedRefinement(model, partition, signatures);
                initializeQuotientData();
                signatures.extendQuotientData(quotientData.value(), options.createQuotientChoiceMapping);
            } else {
                STORM_LOG_THROW_UNCONDITIONALLY(storm::exceptions::NotSupportedException,
                                                "Bisimulation with interval-abstraction and zero tolerance is not supported.");
            }
        } else if constexpr (!std::is_same_v<ValueType, storm::RationalFunction>) {
            storm::bisimulation::Signatures<ValueType, SignatureMode::Approximative, QuotientValueType> signatures(
                model, choiceClasses, partition, storm::utility::convertNumber<storm::IntervalBaseType<ValueType>>(options.tolerance));
            storm::bisimulation::performSignatureBasedRefinement(model, partition, signatures);
            initializeQuotientData();
            signatures.extendQuotientData(quotientData.value(), options.createQuotientChoiceMapping);
        } else {
            STORM_LOG_THROW_UNCONDITIONALLY(storm::exceptions::NotSupportedException, "Bisimulation with positive tolerance "
                                                                                          << storm::utility::convertNumber<double>(options.tolerance)
                                                                                          << " is not supported for parametric models.");
        }
    } else if constexpr (!storm::IsIntervalType<QuotientValueType>) {  // IsIntervalType shouldn't hold either way but condition needed for compilation
        // Use Splitter-based refinement (weak or strong)
        if (isWeak) {
            // Weak bisimulation additionally needs to know the divergent, step sensitive and silent states. Computing them further refines the partition.
            // Both that computation and the refinement need the transposed transition matrix, so we build it once and share it.
            auto const backwardTransitions = model.getBackwardTransitions();
            auto weakData = initialization.getWeakBisimulationData(partition, backwardTransitions, preservationInformation);
            STORM_LOG_STATISTICS(sw << " seconds for weak bisimulation initialization (" << partition.getNumberOfBlocks() << " blocks).");
            sw.restart();
            if (model.isDiscreteTimeModel()) {
                storm::bisimulation::performSplitterBasedRefinement<ValueType, SplitterRefinementMode::WeakDiscreteTime>(
                    model, backwardTransitions, partition, storm::utility::convertNumber<ValueType>(options.tolerance), weakData);
            } else {
                storm::bisimulation::performSplitterBasedRefinement<ValueType, SplitterRefinementMode::WeakContinuousTime>(
                    model, backwardTransitions, partition, storm::utility::convertNumber<ValueType>(options.tolerance), weakData);
            }
            // For weak bisimulation, the quotient transitions have to be derived from a state that can actually leave its block.
            auto const nonSilentStates = ~weakData.silentStates;
            initializeQuotientData(nonSilentStates);
            quotientData->weakData.emplace(std::move(weakData));

        } else {
            storm::bisimulation::performSplitterBasedRefinement<ValueType, SplitterRefinementMode::Strong>(
                model, model.getBackwardTransitions(), partition, storm::utility::convertNumber<ValueType>(options.tolerance));
            initializeQuotientData();
        }
        STORM_LOG_ASSERT(!model.isNondeterministicModel(), "Did not expect nondeterministic model for splitter-based refinement.");
        if (options.createQuotientChoiceMapping) {
            quotientData->toQuotientChoice = quotientData->toQuotientState;  // For deterministic models, quotient state and choice mappings are identical.
        }
    }
    STORM_LOG_ASSERT(quotientData.has_value(), "Quotient data should have been initialized.");  // Assert that initializeQuotientData() has been called.

    // Step 3: Extract the quotient
    auto quotientModel =
        storm::bisimulation::Quotient<ValueType, QuotientValueType>::buildFromPartition(model, options, preservationInformation, quotientData.value());
    STORM_LOG_STATISTICS(sw << " seconds for quotient extraction.");
    STORM_LOG_STATISTICS("-------------------------------------------");

    return {.quotient = std::move(quotientModel),
            .toQuotientStateMapping = std::move(quotientData->toQuotientState),
            .toQuotientChoiceMapping = std::move(quotientData->toQuotientChoice)};
}

// Instantiations with ValueType == QuotientValueType
template ReturnType<double> performBisimulationMinimization(storm::models::sparse::Model<double> const& model,
                                                            std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas, Options const& options);
template ReturnType<storm::RationalNumber> performBisimulationMinimization(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                           std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                           Options const& options);
template ReturnType<storm::RationalFunction> performBisimulationMinimization(storm::models::sparse::Model<storm::RationalFunction> const& model,
                                                                             std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                             Options const& options);
template ReturnType<storm::Interval> performBisimulationMinimization(storm::models::sparse::Model<storm::Interval> const& model,
                                                                     std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                     Options const& options);
template ReturnType<storm::RationalInterval> performBisimulationMinimization(storm::models::sparse::Model<storm::RationalInterval> const& model,
                                                                             std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
                                                                             Options const& options);

// Instantiations for Interval Abstraction: QuotientValueType == IntervalType<ValueType>
template ReturnType<storm::Interval> performBisimulationMinimization<double, storm::Interval>(
    storm::models::sparse::Model<double> const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas, Options const& options);
template ReturnType<storm::RationalInterval> performBisimulationMinimization<storm::RationalNumber, storm::RationalInterval>(
    storm::models::sparse::Model<storm::RationalNumber> const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
    Options const& options);

}  // namespace storm::bisimulation

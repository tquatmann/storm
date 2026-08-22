#include "storm/transformer/bisimulation/Bisimulation.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/models/sparse/Model.h"

#include "storm/transformer/bisimulation/Initialization.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Quotient.h"
#include "storm/transformer/bisimulation/Refinement.h"
#include "storm/transformer/bisimulation/Signatures.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

template<typename ValueType>
ReturnType<ValueType> applyBisimulationMinimization(storm::models::sparse::Model<ValueType> const& model, Options const& options,
                                                    std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) {
    storm::utility::Stopwatch sw(true);
    STORM_LOG_STATISTICS("-------- Bisimulation Minimization --------");
    // Obtain an initial partition based on what needs to be preserved (labels, rewards, ...)
    storm::bisimulation::Initialization<ValueType> initialization(model, options, formulas);
    auto const preservationInformation = initialization.getPreservationInformation();
    auto const choiceClasses = initialization.getChoiceClasses();
    auto partition = initialization.getInitialStatePartition(choiceClasses);
    STORM_LOG_STATISTICS(sw << " seconds for initial partition (" << partition.getNumberOfBlocks() << " blocks).");
    sw.restart();

    // apply refinement using the initial partition, choiceClasses
    if (model.isNondeterministicModel()) {
        if constexpr (!storm::IsIntervalType<ValueType> && !std::is_same_v<ValueType, storm::RationalFunction>) {
            storm::bisimulation::Signatures<ValueType, SignatureMode::Exact> signatures(model, choiceClasses, partition);
            // storm::bisimulation::Signatures<ValueType, SignatureMode::Approximative> signatures(
            //     model, choiceClasses, partition, storm::utility::convertNumber<ValueType>(options.floatTolerance.value_or(0.0)));
            storm::bisimulation::performSignatureBasedRefinement(model, partition, signatures);
            STORM_LOG_STATISTICS(sw << " seconds for signature refinement (" << partition.getNumberOfBlocks() << " blocks).");
        }
    } else {
        ValueType const tolerance = storm::NumberTraits<ValueType>::IsExact ? storm::utility::zero<ValueType>()
                                                                            : storm::utility::convertNumber<ValueType>(options.floatTolerance.value_or(0.0));
        storm::bisimulation::performSplitterBasedRefinement<ValueType>(model, partition, tolerance);
        STORM_LOG_STATISTICS(sw << " seconds for splitter refinement (" << partition.getNumberOfBlocks() << " blocks).");
    }
    sw.restart();

    // extract the quotient
    using ExactQuotient = storm::bisimulation::Quotient<storm::bisimulation::QuotientType::Exact, ValueType>;
    auto const indexMapping = ExactQuotient::computeIndexMappings(model, partition, choiceClasses);
    auto quotientModel = ExactQuotient::buildFromPartition(model, options, preservationInformation, indexMapping);
    STORM_LOG_STATISTICS(sw << " seconds for quotient extraction.");
    STORM_LOG_STATISTICS("-------------------------------------------");

    return {std::move(quotientModel), std::move(indexMapping.toQuotientState)};
}

template ReturnType<double> applyBisimulationMinimization(storm::models::sparse::Model<double> const& model, Options const& options,
                                                          std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);
template ReturnType<storm::RationalNumber> applyBisimulationMinimization(storm::models::sparse::Model<storm::RationalNumber> const& model,
                                                                         Options const& options,
                                                                         std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);
template ReturnType<storm::RationalFunction> applyBisimulationMinimization(storm::models::sparse::Model<storm::RationalFunction> const& model,
                                                                           Options const& options,
                                                                           std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);
template ReturnType<storm::Interval> applyBisimulationMinimization(storm::models::sparse::Model<storm::Interval> const& model, Options const& options,
                                                                   std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);
template ReturnType<storm::RationalInterval> applyBisimulationMinimization(storm::models::sparse::Model<storm::RationalInterval> const& model,
                                                                           Options const& options,
                                                                           std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);

}  // namespace storm::bisimulation

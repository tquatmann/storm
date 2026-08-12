#include "storm/transformer/bisimulation/Bisimulation.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/models/sparse/Model.h"

#include "storm/transformer/bisimulation/Initialization.h"
#include "storm/transformer/bisimulation/Partition.h"
#include "storm/transformer/bisimulation/Quotient.h"
#include "storm/transformer/bisimulation/Refinement.h"
#include "storm/utility/Stopwatch.h"

namespace storm::bisimulation {

template<typename ValueType>
ReturnType<ValueType> applyBisimulationMinimization(storm::models::sparse::Model<ValueType> const& model, Options const& options,
                                                    std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) {
    storm::utility::Stopwatch sw(true);
    // Obtain an initial partition based on what needs to be preserved (labels, rewards, ...)
    storm::bisimulation::Initialization<ValueType> initialization(model, options, formulas);
    auto const preservationInformation = initialization.getPreservationInformation();
    auto const choiceClasses = initialization.getChoiceClasses();
    auto partition = initialization.getInitialStatePartition(choiceClasses);
    STORM_LOG_STATISTICS("Initial partition with " << partition.getNumberOfBlocks() << " blocks computed after " << sw << " seconds.\n");
    sw.restart();

    // apply refinement using the initial partition, choiceClasses
    storm::bisimulation::performPartitionRefinement(model, partition, choiceClasses);
    STORM_LOG_STATISTICS("Refinement terminated with " << partition.getNumberOfBlocks() << " blocks computed after " << sw << " seconds.\n");
    sw.restart();

    // extract the quotient
    using ExactQuotient = storm::bisimulation::Quotient<storm::bisimulation::QuotientType::Exact, ValueType>;
    auto const indexMapping = ExactQuotient::computeIndexMappings(model, partition, choiceClasses);
    auto quotientModel = ExactQuotient::buildFromPartition(model, options, preservationInformation, indexMapping);
    STORM_LOG_STATISTICS("Quotient extracted after " << sw << " seconds.\n");

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

#include <memory>
#include <sstream>
#include <vector>

#include "storm-config.h"
#include "storm-parsers/api/storm-parsers.h"
#include "storm-parsers/parser/AutoParser.h"
#include "storm-parsers/parser/DirectEncodingParser.h"
#include "storm-parsers/parser/FormulaParser.h"
#include "storm/adapters/IntervalAdapter.h"
#include "storm/api/storm.h"
#include "storm/environment/solver/MinMaxSolverEnvironment.h"
#include "storm/modelchecker/prctl/SparseDtmcPrctlModelChecker.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/storage/geometry/Halfspace.h"
#include "storm/storage/stateminimization/bisimulation/DeterministicIntervalModelBisimulationDecomposition.h"
#include "test/storm_gtest.h"

namespace {

std::unique_ptr<storm::modelchecker::QualitativeCheckResult> getInitialStateFilter(
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> const& model) {
    return std::make_unique<storm::modelchecker::ExplicitQualitativeCheckResult<double>>(model->getInitialStates());
}

double getQuantitativeResultAtInitialState(std::shared_ptr<storm::models::sparse::Model<storm::Interval>> const& model,
                                           std::unique_ptr<storm::modelchecker::CheckResult>& result) {
    auto filter = getInitialStateFilter(model);
    result->filter(*filter);
    return result->asQuantitativeCheckResult<double>().getMin();
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, CreatePolytopesfromIDTMC) {
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
        storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/brp-16-2.drn");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(613ul, dtmc->getNumberOfStates());
    EXPECT_TRUE(modelPtr->hasUncertainty());

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeDecomposition());
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(328ul, result->getNumberOfStates());
    EXPECT_EQ(456ul, result->getNumberOfTransitions());

    // Property on BRP: Pmin=? [F (s = 5)]
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, ParseIDTMCFromPrism) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-16-2.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F (s = 5)]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
        storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(bisim.computeDecomposition());
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(328ul, result->getNumberOfStates());
    EXPECT_EQ(456ul, result->getNumberOfTransitions());
}

// TODO: Move these tests into their own class.
// TEST(DeterministicIntervalModelBisimulationDecompositionTest, ApplyEpsilonBisimOnPrism) {
//     std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-32-2.pm";
//     storm::prism::Program program = storm::api::parseProgram(programFile);
//     program = storm::utility::prism::preprocess(program, "");
//     std::string formulasAsString = "Pmin=? [F (s = 5)]";
//     storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
//         storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();
//
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.1);
//
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     EXPECT_EQ(328ul, result->getNumberOfStates());
//     EXPECT_EQ(456ul, result->getNumberOfTransitions());
// }

TEST(DeterministicIntervalModelBisimulationDecompositionTest, Build1_1Interval) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-point-intervals.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "N=16,MAX=2");
    std::string formulasAsString = "Pmin=? [F (s = 5)]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
        storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, IntervalAddition) {
    // Interval.h
    // works using: "using roundingP = boost::numeric::interval_lib::rounded_arith_exact<double>;" or
    // works using: "using roundingP = boost::numeric::interval_lib::rounded_transc_opp<double>;"
    // does not work using (currently standard): "using roundingP =
    // boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double> >;"
    auto first = storm::Interval(0.98, 0.98);
    auto second = storm::Interval(0.02, 0.02);

    auto sum = first + second;

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "sum lower: " << sum.lower() << std::endl;
    std::cout << "sum upper: " << sum.upper() << std::endl;

    EXPECT_TRUE(storm::utility::isOne(sum));
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, IntervalZeroIsEmptySetOnContainerInit) {
    // Initializes elements with (0.0, 0.0) instead of [0.0, 0.0], meaning each interval is the empty set
    std::vector<storm::Interval> cont(1);

    storm::Interval ab(1.0, 2.0);

    // I_1 + I_2 = \{ a + b | a \in I_1, b \in I_2 \}
    auto sum = cont[0] + ab;

    // Fails, as sum is the empty set because cont[0] does not include any elements
    EXPECT_FALSE(sum.lower() == 1.0);

    cont = std::vector<storm::Interval>(1, storm::utility::zero<storm::Interval>());

    sum = cont[0] + ab;

    EXPECT_EQ(sum.lower(), 1.0);
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, IntervalArithmeticSanityCheck) {
    storm::Interval a(0.1, 0.2);
    storm::Interval b(0.0, 0.1);
    storm::Interval c(0.1, 0.2);

    // a + b = [0.1, 0.3]
    EXPECT_DOUBLE_EQ((a + b).lower(), 0.1);
    EXPECT_DOUBLE_EQ((a + b).upper(), 0.3);

    // enhancement: c -> c' = [0.05, 0.2]
    auto enhancedC = (c - storm::Interval(0.05, 0.05)) + storm::Interval(0.0, 0.05);
    EXPECT_DOUBLE_EQ(enhancedC.lower(), 0.05);
    EXPECT_DOUBLE_EQ(enhancedC.upper(), 0.2);
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny02IDTMC) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/tiny-02.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F \"a\"]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
        storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeDecomposition());
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(2ul, result->getNumberOfStates());
    EXPECT_EQ(2ul, result->getNumberOfTransitions());

    result->printModelInformationToStream(std::cout);

    // Property on BRP: Pmin=? [F (s = 5)]
}

// TODO: Create quotient of toy example manually on paper and compare
TEST(DeterministicIntervalModelBisimulationDecompositionTest, TinyIDTMC) {
    // TODO: Parser has a bug. It does not check whether the drn-file represents a valid IDTMC.
    // TODO: Set transition: 4 : [0.1, 0.3] of state 4 to [0.1, 0.2], then there is no possible realization s.t. p_1 + p_2 = 1

    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
        storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/tiny-05.drn");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(7ul, dtmc->getNumberOfStates());
    EXPECT_TRUE(modelPtr->hasUncertainty());

    modelPtr->printModelInformationToStream(std::cout);

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeDecomposition());
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(4ul, result->getNumberOfStates());
    EXPECT_EQ(6ul, result->getNumberOfTransitions());

    result->printModelInformationToStream(std::cout);

    // Property on BRP: Pmin=? [F (s = 5)]
}

// TODO: Move these tests into their own class.
// TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny03IDTMC) {
//     std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/tiny-03.pm";
//     storm::prism::Program program = storm::api::parseProgram(programFile);
//     program = storm::utility::prism::preprocess(program, "");
//     std::string formulasAsString = "Pmin=? [F (s = 3)];Pmin=? [F (s = 1)]";
//     std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
//         storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
//         storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();
//
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.8);
//
//     dtmc->printModelInformationToStream(std::cout);
//
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     EXPECT_EQ(3ul, result->getNumberOfStates());
//     EXPECT_EQ(5ul, result->getNumberOfTransitions());
// }

// TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny04IDTMC_Paper) {
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
//         storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/tiny-04-paper.drn");
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
//     ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
//     ASSERT_EQ(4ul, dtmc->getNumberOfStates());
//     ASSERT_EQ(9ul, dtmc->getNumberOfTransitions());
//     EXPECT_TRUE(modelPtr->hasUncertainty());
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.080001);
//     modelPtr->printModelInformationToStream(std::cout);
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     EXPECT_EQ(2ul, result->getNumberOfStates());
//     EXPECT_EQ(3ul, result->getNumberOfTransitions());
//     result->printModelInformationToStream(std::cout);
// }

// TEST(DeterministicIntervalModelBisimulationDecompositionTest, CrowdsIDTMC_Perturbed_RationalInterval) {
//     std::shared_ptr<storm::models::sparse::Model<storm::RationalInterval>> modelPtr =
//         storm::parser::parseDirectEncodingModel<storm::RationalInterval>(STORM_TEST_RESOURCES_DIR "/idtmc/after_perturbation_0.000001_2026-05-05_16-11-37");
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalInterval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::RationalInterval>>();
//     ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
//     EXPECT_TRUE(modelPtr->hasUncertainty());
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>>::BisimulationOptions
//         options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.001);
//     modelPtr->printModelInformationToStream(std::cout);
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::RationalInterval>> bisim(*dtmc, options);
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     std::shared_ptr<storm::models::sparse::Model<storm::RationalInterval>> result;
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     result->printModelInformationToStream(std::cout);
// }

// TEST(DeterministicIntervalModelBisimulationDecompositionTest, CrowdsIDTMC_Perturbed_DoubleInterval_Epsilon) {
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
//         storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/after_perturbation_0.000001_2026-05-05_16-11-37");
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
//     ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
//     EXPECT_TRUE(modelPtr->hasUncertainty());
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.001);
//     modelPtr->printModelInformationToStream(std::cout);
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     result->printModelInformationToStream(std::cout);
// }

TEST(DeterministicIntervalModelBisimulationDecompositionTest, CrowdsIDTMC_Perturbed_DoubleInterval_Exact) {
    storm::Interval intervalA = storm::Interval(1, carl::BoundType::WEAK, 1, carl::BoundType::WEAK);
    storm::Interval intervalB = storm::Interval(storm::utility::zero<storm::IntervalBaseType<storm::Interval>>(), carl::BoundType::WEAK,
                                                storm::utility::one<storm::IntervalBaseType<storm::Interval>>(), carl::BoundType::WEAK);

    // TODO: check exact bisimulation result for IDTMCs. compare with non-interval model. check quotient for weird intervals.
    // TODO: try to get --exact to run. maybe merge master into this branch again?

    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
        storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/after_perturbation_0.000001_2026-05-05_16-11-37");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    EXPECT_TRUE(modelPtr->hasUncertainty());
    typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
    // "P=? [F \"observe0Greater1\" ]"
    options.preserveFormula(*storm::api::extractFormulasFromProperties(storm::api::parseProperties("P=? [F \"observe0Greater1\" ]")).at(0).get());
    modelPtr->printModelInformationToStream(std::cout);
    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
    ASSERT_NO_THROW(bisim.computeDecomposition());
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(result = bisim.getQuotient());
    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    result->printModelInformationToStream(std::cout);
}

// TODO: Move these tests into their own class.
// TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny04IDTMC_Paper_PM_Epsilon) {
//     std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/tiny-04-paper.pm";
//     storm::prism::Program program = storm::api::parseProgram(programFile);
//     program = storm::utility::prism::preprocess(program, "");
//     std::string formulasAsString = "Pmin=? [F (s = 3)]";
//     std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
//         storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
//
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
//         storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();
//
//     ASSERT_EQ(storm::models::ModelType::Dtmc, dtmc->getType());
//     ASSERT_EQ(4ul, dtmc->getNumberOfStates());
//     ASSERT_EQ(9ul, dtmc->getNumberOfTransitions());
//     EXPECT_TRUE(dtmc->hasUncertainty());
//
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.080001);  // epsilon 0.08 should make this test work, floating point error requires slightly more
//
//     dtmc->printModelInformationToStream(std::cout);
//
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//
//     result->printModelInformationToStream(std::cout);
//
//     storm::api::exportSparseModelAsDrn(result, "after_bisimulation_test", std::vector<std::string>(), false);
//
//     EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
//     EXPECT_EQ(2ul, result->getNumberOfStates());
//     EXPECT_EQ(3ul, result->getNumberOfTransitions());
// }

TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny04IDTMC_Paper_PM) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/tiny-04-paper.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F (s = 3)]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));

    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
        storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();

    ASSERT_EQ(storm::models::ModelType::Dtmc, dtmc->getType());
    ASSERT_EQ(4ul, dtmc->getNumberOfStates());
    ASSERT_EQ(9ul, dtmc->getNumberOfTransitions());
    EXPECT_TRUE(dtmc->hasUncertainty());

    typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
    // options.preserveFormula(*formulas[0].get());

    dtmc->printModelInformationToStream(std::cout);

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(bisim.computeDecomposition());
    ASSERT_NO_THROW(result = bisim.getQuotient());

    result->printModelInformationToStream(std::cout);

    storm::api::exportSparseModelAsDrn(result, "after_bisimulation_test", std::vector<std::string>(), false);

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(3ul, result->getNumberOfStates());
    EXPECT_EQ(6ul, result->getNumberOfTransitions());
}

// TODO: Move these tests into their own class.
// TEST(DeterministicIntervalModelBisimulationDecompositionTest, BRP_Quotient_MC) {
//     std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-point-intervals.pm";
//     storm::prism::Program program = storm::api::parseProgram(programFile);
//     program = storm::utility::prism::preprocess(program, "N=16,MAX=2");
//     std::string formulasAsString = "P=? [ F s=5 ]";
//     std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
//         storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc =
//         storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Dtmc<storm::Interval>>();
//
//     typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>::BisimulationOptions options;
//     // options.preserveFormula(*formulas[0].get());
//     options.setUsesEpsilon(true);
//     options.setEpsilon(0.02);  // epsilon 0.08 should make this test work, floating point error requires slightly more
//
//     dtmc->printModelInformationToStream(std::cout);
//
//     storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc, options);
//     std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
//     ASSERT_NO_THROW(bisim.computeDecomposition());
//     ASSERT_NO_THROW(result = bisim.getQuotient());
//
//     storm::Environment env;
//     env.solver().minMax().setMethod(storm::solver::MinMaxMethod::ValueIteration);
//
//     std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> quotient = result->as<storm::models::sparse::Dtmc<storm::Interval>>();
//     ASSERT_EQ(storm::models::ModelType::Dtmc, quotient->getType());
//     auto taskMin = storm::modelchecker::CheckTask<storm::logic::Formula, double>(*formulas[0]);
//     taskMin.setUncertaintyResolutionMode(storm::UncertaintyResolutionMode::Minimize);
//
//     auto checker = storm::modelchecker::SparseDtmcPrctlModelChecker<storm::models::sparse::Dtmc<storm::Interval>>(*quotient);
//     auto resultMin = checker.check(env, taskMin);
//
//     std::cout << getQuantitativeResultAtInitialState(quotient, resultMin) << std::endl;
// }

}  // namespace
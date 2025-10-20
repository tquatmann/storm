#include <sstream>
#include "storm/models/sparse/Dtmc.h"

#include "test/storm_gtest.h"

#include <storm-parsers/parser/DirectEncodingParser.h>
#include <memory>
#include <storm/storage/bisimulation/DeterministicIntervalModelBisimulationDecomposition.h>
#include <vector>
#include "storm/storage/geometry/Halfspace.h"

#include "storm/modelchecker/prctl/SparseDtmcPrctlModelChecker.h"

#include "storm-config.h"
#include "storm-parsers/parser/AutoParser.h"
#include "storm-parsers/parser/FormulaParser.h"
#include "storm/api/storm.h"

#include "storm-parsers/api/storm-parsers.h"
#include "storm/api/storm.h"

namespace {

TEST(DeterministicIntervalModelBisimulationDecompositionTest, CreatePolytopesfromIDTMC) {
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
        storm::parser::DirectEncodingParser<storm::Interval>::parseModel(STORM_TEST_RESOURCES_DIR "/idtmc/brp-16-2.drn");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(613ul, dtmc->getNumberOfStates());
    EXPECT_TRUE(modelPtr->hasUncertainty());

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeBisimulationDecomposition());
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
    std::shared_ptr<storm::models::sparse::Dtmc<carl::Interval<double>>> dtmc =
        storm::api::buildSparseModel<carl::Interval<double>>(program, formulas)->as<storm::models::sparse::Dtmc<carl::Interval<double>>>();

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<carl::Interval<double>>> bisim(*dtmc);
    std::shared_ptr<storm::models::sparse::Model<carl::Interval<double>>> result;
    ASSERT_NO_THROW(bisim.computeBisimulationDecomposition());
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(328ul, result->getNumberOfStates());
    EXPECT_EQ(456ul, result->getNumberOfTransitions());
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, ApplyEpsilonBisimOnPrism) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-32-2.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F (s = 5)]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<carl::Interval<double>>> dtmc =
        storm::api::buildSparseModel<carl::Interval<double>>(program, formulas)->as<storm::models::sparse::Dtmc<carl::Interval<double>>>();

    typename storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<carl::Interval<double>>>::Options options;
    options.preserveFormula(*formulas[0].get());
    options.setUsesEpsilon(true);
    options.setEpsilon(0.1);

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<carl::Interval<double>>> bisim(*dtmc, options);
    std::shared_ptr<storm::models::sparse::Model<carl::Interval<double>>> result;
    ASSERT_NO_THROW(bisim.computeBisimulationDecomposition());
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(328ul, result->getNumberOfStates());
    EXPECT_EQ(456ul, result->getNumberOfTransitions());
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, Build1_1Interval) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/brp-32-2-point-intervals.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F (s = 5)]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<carl::Interval<double>>> dtmc =
        storm::api::buildSparseModel<carl::Interval<double>>(program, formulas)->as<storm::models::sparse::Dtmc<carl::Interval<double>>>();
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, IntervalAddition) {
    auto first = storm::Interval(0.98, 0.98);
    auto second = storm::Interval(0.02, 0.02);

    auto sum = first + second;

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "sum lower: " << sum.lower() << std::endl;
    std::cout << "sum upper: " << sum.upper() << std::endl;

    EXPECT_TRUE(storm::utility::isOne(sum));
}

TEST(DeterministicIntervalModelBisimulationDecompositionTest, Tiny02IDTMC) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/idtmc/tiny-02.pm";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmin=? [F \"a\"]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Dtmc<carl::Interval<double>>> dtmc =
        storm::api::buildSparseModel<carl::Interval<double>>(program, formulas)->as<storm::models::sparse::Dtmc<carl::Interval<double>>>();

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeBisimulationDecomposition());
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
        storm::parser::DirectEncodingParser<storm::Interval>::parseModel(STORM_TEST_RESOURCES_DIR "/idtmc/tiny-01.drn");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(7ul, dtmc->getNumberOfStates());
    EXPECT_TRUE(modelPtr->hasUncertainty());

    modelPtr->printModelInformationToStream(std::cout);

    storm::storage::DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>> bisim(*dtmc);
    ASSERT_NO_THROW(bisim.computeBisimulationDecomposition());
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(result = bisim.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Dtmc, result->getType());
    EXPECT_EQ(4ul, result->getNumberOfStates());
    EXPECT_EQ(6ul, result->getNumberOfTransitions());

    result->printModelInformationToStream(std::cout);

    // Property on BRP: Pmin=? [F (s = 5)]
}

}
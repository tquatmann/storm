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

// TODO: Create quotient of toy example manually on paper and compare
TEST(DeterministicIntervalModelBisimulationDecompositionTest, TinyIDTMC) {
    // TODO: Parser has a bug. It does not check whether the drn-file represents a valid IDTMC.
    // TODO: Set transition: 4 : [0.1, 0.3] of state 4 to [0.1, 0.2], then there is no possible realization s.t. p_1 + _2 = 1

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
#include <sstream>
#include "storm/models/sparse/Dtmc.h"

#include "test/storm_gtest.h"

#include <storm-parsers/parser/DirectEncodingParser.h>
#include <memory>
#include <storm/storage/bisimulation/DeterministicIntervalModelBisimulationDecomposition.h>
#include <vector>
#include "storm/storage/geometry/Halfspace.h"

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
    ASSERT_NO_THROW(bisim.getQuotient());
}

// TODO: Create test case for toy example on Remarkable

}
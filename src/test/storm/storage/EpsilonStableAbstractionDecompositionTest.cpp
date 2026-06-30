#include <memory>
#include <sstream>
#include <vector>

#include "storm-config.h"
#include "storm-parsers/api/storm-parsers.h"
#include "storm-parsers/parser/AutoParser.h"
#include "storm-parsers/parser/FormulaParser.h"
#include "storm/adapters/IntervalAdapter.h"
#include "storm/api/storm.h"
#include "storm/storage/stateminimization/bisimulation/IntervalModelBisimulationDecomposition.h"
#include "storm/storage/umb/import/SparseModelFromUmb.h"
#include "storm/storage/umb/import/UmbImport.h"
#include "storm/storage/umb/model/UmbModel.h"
#include "test/storm_gtest.h"

namespace {
TEST(EpsilonStableAbstractionDecompositionTest, AircraftTinyImdpPointIntervals) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/imdp/aircraft-tiny-point-intervals.prism";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmax=? [!collision U \"goal\"]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Mdp<storm::Interval>> mdp =
        storm::api::buildSparseModel<storm::Interval>(program, formulas)->as<storm::models::sparse::Mdp<storm::Interval>>();

    typename storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>>::EpsilonStableAbstractionOptions
        options;
    options.preserveFormula(*formulas[0].get());
    options.setEpsilon(0.001);

    storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>> abstractionDecomposition(*mdp, options);
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(abstractionDecomposition.computeDecomposition());
    ASSERT_NO_THROW(result = abstractionDecomposition.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Mdp, result->getType());
    EXPECT_EQ(86ul, result->getNumberOfStates());
    EXPECT_EQ(798ul, result->getNumberOfTransitions());
    EXPECT_EQ(254, result->getNumberOfChoices());
}

TEST(EpsilonStableAbstractionDecompositionTest, UAV2DUmbImdpIntervals) {
    storm::umb::ImportOptions importOptions;
    importOptions.buildStateValuations = false;
    auto umb = storm::umb::importUmb(std::filesystem::path{STORM_TEST_RESOURCES_DIR "/imdp/uav-2d.umb"}, importOptions);

    std::shared_ptr<storm::models::sparse::Mdp<storm::Interval>> model =
        storm::umb::sparseModelFromUmb<storm::Interval>(umb, importOptions)->as<storm::models::sparse::Mdp<storm::Interval>>();

    typename storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>>::EpsilonStableAbstractionOptions
        options;
    options.setEpsilon(0.001);

    storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>> abstractionDecomposition(*model, options);
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(abstractionDecomposition.computeDecomposition());
    ASSERT_NO_THROW(result = abstractionDecomposition.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Mdp, result->getType());
    EXPECT_EQ(190ul, result->getNumberOfStates());
    EXPECT_EQ(78118ul, result->getNumberOfTransitions());
    EXPECT_EQ(3679ul, result->getNumberOfChoices());
}

TEST(EpsilonStableAbstractionDecompositionTest, TestLearning) {
    std::string programFile = STORM_TEST_RESOURCES_DIR "/mdp/aircraft-tiny.prism";
    storm::prism::Program program = storm::api::parseProgram(programFile);
    program = storm::utility::prism::preprocess(program, "");
    std::string formulasAsString = "Pmax=? [!collision U \"goal\"]";
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas =
        storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
    std::shared_ptr<storm::models::sparse::Mdp<double>> mdp = storm::api::buildSparseModel<double>(program, formulas)->as<storm::models::sparse::Mdp<double>>();

    storm::api::learnIMDPFromMDPByClopperPearson<double>(mdp, 0.01, 10);
}

TEST(EpsilonStableAbstractionDecompositionTest, TestFirewire) {
    storm::umb::ImportOptions importOptions;
    importOptions.buildStateValuations = false;
    auto umb = storm::umb::importUmb(std::filesystem::path{STORM_TEST_RESOURCES_DIR "/imdp/firewire_exact_point.umb"}, importOptions);

    std::shared_ptr<storm::models::sparse::Mdp<storm::Interval>> model =
        storm::umb::sparseModelFromUmb<storm::Interval>(umb, importOptions)->as<storm::models::sparse::Mdp<storm::Interval>>();

    typename storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>>::EpsilonStableAbstractionOptions
        options;
    options.setEpsilon(0.001);

    storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<storm::Interval>> abstractionDecomposition(*model, options);
    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> result;
    ASSERT_NO_THROW(abstractionDecomposition.computeDecomposition());
    ASSERT_NO_THROW(result = abstractionDecomposition.getQuotient());

    EXPECT_EQ(storm::models::ModelType::Mdp, result->getType());
    EXPECT_EQ(190ul, result->getNumberOfStates());
    EXPECT_EQ(78118ul, result->getNumberOfTransitions());
    EXPECT_EQ(3679ul, result->getNumberOfChoices());
}
}  // namespace

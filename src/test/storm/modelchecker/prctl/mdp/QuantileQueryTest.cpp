#include "storm-config.h"
#include "test/storm_gtest.h"

#include "test/storm_gtest.h"

#include "storm-parsers/api/model_descriptions.h"
#include "storm-parsers/api/properties.h"
#include "storm/api/builder.h"
#include "storm/api/properties.h"
#include "storm/parser/CSVParser.h"

#include "storm/environment/solver/MinMaxSolverEnvironment.h"
#include "storm/logic/Formulas.h"
#include "storm/modelchecker/prctl/SparseDtmcPrctlModelChecker.h"
#include "storm/modelchecker/prctl/SparseMdpPrctlModelChecker.h"
#include "storm/modelchecker/results/ExplicitParetoCurveCheckResult.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/jani/Property.h"

namespace {

enum class MdpEngine { PrismSparse, JaniSparse, Hybrid, PrismDd, JaniDd };

class UnsoundEnvironment {
   public:
    typedef double ValueType;
    static storm::Environment createEnvironment() {
        storm::Environment env;
        env.solver().minMax().setPrecision(storm::utility::convertNumber<storm::RationalNumber>(1e-10));
        return env;
    }
};

class SoundEnvironment {
   public:
    typedef double ValueType;
    static storm::Environment createEnvironment() {
        storm::Environment env;
        env.solver().setForceSoundness(true);
        return env;
    }
};

class ExactEnvironment {
   public:
    typedef storm::RationalNumber ValueType;
    static storm::Environment createEnvironment() {
        storm::Environment env;
        return env;
    }
};

template<typename TestType>
class QuantileQueryTest : public ::testing::Test {
   public:
    typedef typename TestType::ValueType ValueType;
    QuantileQueryTest() : _environment(TestType::createEnvironment()) {}
    storm::Environment const& env() const {
        return _environment;
    }
    ValueType parseNumber(std::string const& input) const {
        if (input.find("inf") != std::string::npos) {
            return storm::utility::infinity<ValueType>();
        }
        return storm::utility::convertNumber<ValueType>(input);
    }

    template<typename MT>
    std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>> buildModelFormulas(
        std::string const& pathToPrismFile, std::string const& formulasAsString, std::string const& constantDefinitionString = "") const {
        std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>> result;
        storm::prism::Program program = storm::api::parseProgram(pathToPrismFile);
        program = storm::utility::prism::preprocess(program, constantDefinitionString);
        result.second = storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
        result.first = storm::api::buildSparseModel<ValueType>(program, result.second)->template as<MT>();
        return result;
    }

    std::vector<storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>> getTasks(
        std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) const {
        std::vector<storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>> result;
        for (auto const& f : formulas) {
            result.emplace_back(*f, true);  // only initial states are relevant
        }
        return result;
    }

    template<typename MT>
    typename std::enable_if<std::is_same<MT, storm::models::sparse::Mdp<ValueType>>::value,
                            std::shared_ptr<storm::modelchecker::AbstractModelChecker<MT>>>::type
    createModelChecker(std::shared_ptr<MT> const& model) const {
        return std::make_shared<storm::modelchecker::SparseMdpPrctlModelChecker<MT>>(*model);
    }
    template<typename MT>
    typename std::enable_if<std::is_same<MT, storm::models::sparse::Dtmc<ValueType>>::value,
                            std::shared_ptr<storm::modelchecker::AbstractModelChecker<MT>>>::type
    createModelChecker(std::shared_ptr<MT> const& model) const {
        return std::make_shared<storm::modelchecker::SparseDtmcPrctlModelChecker<MT>>(*model);
    }

    std::pair<bool, std::string> compareResult(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model,
                                               std::unique_ptr<storm::modelchecker::CheckResult>& result, std::vector<std::string> const& expected) {
        bool equal = true;
        std::string errorMessage = "";
        ValueType comparePrecision =
            std::is_same<ValueType, double>::value ? storm::utility::convertNumber<ValueType>(1e-10) : storm::utility::zero<ValueType>();
        auto filter = getInitialStateFilter(model);
        result->filter(*filter);
        std::vector<std::vector<ValueType>> resultPoints;
        if (result->isExplicitParetoCurveCheckResult()) {
            resultPoints = result->asExplicitParetoCurveCheckResult<ValueType>().getPoints();
        } else {
            if (!result->isExplicitQuantitativeCheckResult()) {
                return {false, "The CheckResult has unexpected type."};
            }
            resultPoints = {{result->asExplicitQuantitativeCheckResult<ValueType>().getMax()}};
        }
        std::vector<std::vector<ValueType>> expectedPoints;
        for (auto const& pointAsString : expected) {
            std::vector<ValueType> point;
            for (auto const& entry : storm::parser::parseCommaSeperatedValues(pointAsString)) {
                point.push_back(parseNumber(entry));
            }
            expectedPoints.push_back(std::move(point));
        }
        for (auto const& resPoint : resultPoints) {
            bool contained = false;
            for (auto const& expPoint : expectedPoints) {
                if (storm::utility::vector::equalModuloPrecision(resPoint, expPoint, comparePrecision, true)) {
                    contained = true;
                    break;
                }
            }
            if (!contained) {
                equal = false;
                errorMessage += "The result " + storm::utility::vector::toString(resPoint) + " is not contained in the expected solution. ";
            }
        }
        for (auto const& expPoint : expectedPoints) {
            bool contained = false;
            for (auto const& resPoint : resultPoints) {
                if (storm::utility::vector::equalModuloPrecision(resPoint, expPoint, comparePrecision, true)) {
                    contained = true;
                    break;
                }
            }
            if (!contained) {
                equal = false;
                errorMessage += "The expected " + storm::utility::vector::toString(expPoint) + " is not contained in the obtained solution. ";
            }
        }
        return {equal, errorMessage};
    }

   private:
    storm::Environment _environment;

    std::unique_ptr<storm::modelchecker::QualitativeCheckResult> getInitialStateFilter(
        std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model) const {
        return std::make_unique<storm::modelchecker::ExplicitQualitativeCheckResult>(model->getInitialStates());
    }
};

typedef ::testing::Types<UnsoundEnvironment, SoundEnvironment, ExactEnvironment> TestingTypes;

TYPED_TEST_SUITE(QuantileQueryTest, TestingTypes, );

TYPED_TEST(QuantileQueryTest, simple_Dtmc) {
    typedef storm::models::sparse::Dtmc<typename TestFixture::ValueType> ModelType;

    std::string formulasString = "quantile(max A, max B, P>0.95 [F{\"first\"}<=A,{\"second\"}<=B s=3]);\n";

    auto modelFormulas = this->template buildModelFormulas<ModelType>(STORM_TEST_RESOURCES_DIR "/dtmc/quantiles_simple_dtmc.pm", formulasString);
    auto model = std::move(modelFormulas.first);
    auto tasks = this->getTasks(modelFormulas.second);
    EXPECT_EQ(4ul, model->getNumberOfStates());
    EXPECT_EQ(6ul, model->getNumberOfTransitions());
    auto checker = this->template createModelChecker<ModelType>(model);
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    std::vector<std::string> expectedResult;
    std::pair<bool, std::string> compare;
    uint64_t taskId = 0;

    expectedResult.clear();
    expectedResult.push_back("7, 18");
    expectedResult.push_back("8, 16");
    expectedResult.push_back("9, 14");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;
}

TYPED_TEST(QuantileQueryTest, simple_Mdp) {
    typedef storm::models::sparse::Mdp<typename TestFixture::ValueType> ModelType;

    std::string formulasString = "quantile(B1, B2, Pmax>0.5 [F{\"first\"}>=B1,{\"second\"}>=B2 s=1]);\n";
    formulasString += "quantile(B1, B2, Pmax>0.5 [F{\"first\"}<=B1,{\"second\"}<=B2 s=1]);\n";
    formulasString += "quantile(B1, B2, Pmin>0.5 [F{\"first\"}<=B1,{\"second\"}<=B2 s=1]);\n";
    formulasString += "quantile(B1, Pmax>0.5 [F{\"second\"}>=B1 s=1]);\n";
    formulasString += "quantile(B1, B2, B3, Pmax>0.5 [F{\"first\"}<=B1,{\"second\"}<=B2,{\"third\"}<=B3 s=1]);\n";

    auto modelFormulas = this->template buildModelFormulas<ModelType>(STORM_TEST_RESOURCES_DIR "/mdp/quantiles_simple_mdp.nm", formulasString);
    auto model = std::move(modelFormulas.first);
    auto tasks = this->getTasks(modelFormulas.second);
    EXPECT_EQ(7ul, model->getNumberOfStates());
    EXPECT_EQ(13ul, model->getNumberOfTransitions());
    auto checker = this->template createModelChecker<ModelType>(model);
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    std::vector<std::string> expectedResult;
    std::pair<bool, std::string> compare;
    uint64_t taskId = 0;

    expectedResult.clear();
    expectedResult.push_back("  0, 6");
    expectedResult.push_back("  1, 5");
    expectedResult.push_back("  2, 4");
    expectedResult.push_back("  3, 3");
    expectedResult.push_back("inf, 2");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;

    expectedResult.clear();
    expectedResult.push_back("  0, 2");
    expectedResult.push_back("  5, 1");
    expectedResult.push_back("  6, 0");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;

    expectedResult.clear();
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;

    expectedResult.clear();
    expectedResult.push_back("  6");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;

    expectedResult.clear();
    expectedResult.push_back("  0, 2, 1.4");
    expectedResult.push_back("  5, 1, 9.8");
    expectedResult.push_back("  6, 0, 9.8");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;
}

TYPED_TEST(QuantileQueryTest, firewire) {
    typedef storm::models::sparse::Mdp<typename TestFixture::ValueType> ModelType;

    std::string formulasString = "quantile(min TIME, min ROUNDS, Pmax>0.95 [F{\"time\"}<=TIME,{\"rounds\"}<=ROUNDS \"done\"]);\n";

    auto modelFormulas = this->template buildModelFormulas<ModelType>(STORM_TEST_RESOURCES_DIR "/mdp/quantiles_firewire.nm", formulasString, "delay=36");
    auto model = std::move(modelFormulas.first);
    auto tasks = this->getTasks(modelFormulas.second);
    EXPECT_EQ(776ul, model->getNumberOfStates());
    EXPECT_EQ(1411ul, model->getNumberOfTransitions());
    auto checker = this->template createModelChecker<ModelType>(model);
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    std::vector<std::string> expectedResult;
    std::pair<bool, std::string> compare;
    uint64_t taskId = 0;

    expectedResult.clear();
    expectedResult.push_back("123, 1");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;
}

TYPED_TEST(QuantileQueryTest, resources) {
    typedef storm::models::sparse::Mdp<typename TestFixture::ValueType> ModelType;

    std::string formulasString = "quantile(max GOLD, max GEM, Pmax>0.95 [F{\"gold\"}>=GOLD,{\"gem\"}>=GEM,{\"steps\"}<=100 true]);\n";

    auto modelFormulas = this->template buildModelFormulas<ModelType>(STORM_TEST_RESOURCES_DIR "/mdp/quantiles_resources.nm", formulasString);
    auto model = std::move(modelFormulas.first);
    auto tasks = this->getTasks(modelFormulas.second);
    EXPECT_EQ(94ul, model->getNumberOfStates());
    EXPECT_EQ(326ul, model->getNumberOfTransitions());
    auto checker = this->template createModelChecker<ModelType>(model);
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    std::vector<std::string> expectedResult;
    std::pair<bool, std::string> compare;
    uint64_t taskId = 0;

    expectedResult.clear();
    expectedResult.push_back("0, 10");
    expectedResult.push_back("1, 9");
    expectedResult.push_back("4, 8");
    expectedResult.push_back("7, 7");
    expectedResult.push_back("8, 4");
    expectedResult.push_back("9, 2");
    expectedResult.push_back("10, 0");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;
}

template<typename ValueType>
auto buildTransitionRewardsTestModel() {
    // There currently is no nice way to input models with transition branch rewards. We're thus doing this manually here.
    storm::storage::SparseMatrixBuilder<ValueType> transitionsBuilder(6, 5, 10, true, true, 5), rewards1Builder(6, 5, 4, true, true, 5),
        rewards2Builder(6, 5, 2, true, true, 5);
    uint64_t row = 0;
    auto newRowGroup = [&]() {
        transitionsBuilder.newRowGroup(row);
        rewards1Builder.newRowGroup(row);
        rewards2Builder.newRowGroup(row);
    };
    auto addNextValue = [&](uint64_t column, std::string probability, std::string reward1, std::string reward2) {
        transitionsBuilder.addNextValue(row, column, storm::utility::convertNumber<ValueType>(probability));
        if (!reward1.empty()) {
            rewards1Builder.addNextValue(row, column, storm::utility::convertNumber<ValueType>(reward1));
        }
        if (!reward2.empty()) {
            rewards2Builder.addNextValue(row, column, storm::utility::convertNumber<ValueType>(reward2));
        }
    };
    newRowGroup();
    addNextValue(0, "99/100", "1", {});
    addNextValue(1, "1/100", "1", "10");
    ++row;
    addNextValue(1, "4/10", "1", "20");
    addNextValue(2, "6/10", {}, {});
    ++row;
    newRowGroup();
    addNextValue(3, "1/2", {}, {});
    addNextValue(4, "1/2", {}, {});
    ++row;
    newRowGroup();
    addNextValue(0, "3/10", {}, {});
    addNextValue(1, "7/10", "3", {});
    ++row;
    newRowGroup();
    addNextValue(3, "1", {}, {});
    ++row;
    newRowGroup();
    addNextValue(4, "1", {}, {});

    storm::models::sparse::StateLabeling labels(5);
    labels.addLabel("init");
    labels.addLabelToState("init", 0);
    labels.addLabel("target");
    labels.addLabelToState("target", 3);

    auto mdp = std::make_shared<storm::models::sparse::Mdp<ValueType>>(transitionsBuilder.build(), std::move(labels));

    std::vector<ValueType> actionRewards1(6);
    actionRewards1[2] = storm::utility::convertNumber<ValueType, uint64_t>(30);
    actionRewards1[3] = storm::utility::convertNumber<ValueType, uint64_t>(2);
    mdp->addRewardModel("rew1", storm::models::sparse::StandardRewardModel<ValueType>(std::nullopt, std::move(actionRewards1), rewards1Builder.build()));
    mdp->addRewardModel("rew2", storm::models::sparse::StandardRewardModel<ValueType>(std::nullopt, std::nullopt, rewards2Builder.build()));
    return mdp;
    // Equivalent PRISM formulation with intermediate states
    // mdp
    //
    // module main
    //
    // x : [0..4];
    // xNext : [0..4];
    //
    //     [tick1] x=0 -> 0.99: (xNext'=0) + 0.01: (xNext'=1);
    //     [tick2] x=0 -> 0.4: (xNext'=1) + 0.6: (xNext'=2);
    //     [tick1] x=1 -> 0.5: (xNext'=3) + 0.5: (xNext'=4);
    //     [tick1] x=2 -> 0.3: (xNext'=0) + 0.7: (xNext'=1);
    //
    //[tack] true -> 1: (x'=xNext);
    // endmodule
    //
    // module ticktack
    //
    // act : [0..2] init 0;
    //
    //[tick1] act=0 -> 1: (act'=1);
    //[tick2] act=0 -> 1: (act'=2);
    //[tack] act>0 -> 1: (act'=0);
    //
    // endmodule
    //
    // rewards "rew2"
    //[tack] x=0 & xNext=1 & act=1: 10;
    //[tack] x=0 & xNext=1 & act=2: 20;
    // endrewards
    //
    // rewards "rew1"
    //[tack] x=0 & act=1: 1;
    //[tack] x=0 & xNext=1 & act=2: 1;
    //[tack] x=1: 30;
    //[tack] x=2: (xNext=1 ? 5 : 2);
    // endrewards
    //
    // label "target" = x=3;
}

TYPED_TEST(QuantileQueryTest, transition_rewards) {
    using ValueType = typename TestFixture::ValueType;
    typedef storm::models::sparse::Mdp<ValueType> ModelType;
    std::string formulasString = "quantile(min B2, Pmax>0.49 [F{\"rew2\"}<=B2 \"target\"]);\n";
    formulasString += "quantile(min B1, min B2, Pmax>0.4 [F{\"rew1\"}<=B1,{\"rew2\"}<=B2 \"target\"]);\n";
    auto formulas = storm::api::extractFormulasFromProperties(storm::api::parseProperties(formulasString));
    auto tasks = this->getTasks(formulas);

    auto model = buildTransitionRewardsTestModel<ValueType>();
    auto checker = this->template createModelChecker<ModelType>(model);
    std::unique_ptr<storm::modelchecker::CheckResult> result;
    std::vector<std::string> expectedResult;
    std::pair<bool, std::string> compare;
    uint64_t taskId = 0;

    expectedResult.clear();
    expectedResult.push_back("10");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;

    expectedResult.clear();
    expectedResult.push_back("129, 10");
    expectedResult.push_back("35, 20");
    result = checker->check(this->env(), tasks[taskId++]);
    compare = this->compareResult(model, result, expectedResult);
    EXPECT_TRUE(compare.first) << compare.second;
}

}  // namespace

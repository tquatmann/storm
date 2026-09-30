#include "BisimulationTestHelper.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/exceptions/NotSupportedException.h"

namespace {

using storm::test::bisimulation::buildModel;
using storm::test::bisimulation::Options;
using storm::test::bisimulation::strongOptions;
using storm::test::bisimulation::weakOptions;

/*!
 * The states 0 and 1 move to the states 1 and 2 with the same (interval) probabilities, so they are bisimilar. State 2 is absorbing and labeled "goal".
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Dtmc<ValueType>> buildDtmc(ValueType const& value) {
    storm::storage::SparseMatrixBuilder<ValueType> builder(3, 3);
    builder.addNextValue(0, 1, value);
    builder.addNextValue(0, 2, storm::utility::one<ValueType>() - value);
    builder.addNextValue(1, 1, value);
    builder.addNextValue(1, 2, storm::utility::one<ValueType>() - value);
    builder.addNextValue(2, 2, storm::utility::one<ValueType>());
    return buildModel<storm::models::sparse::Dtmc<ValueType>>(builder.build(), {{"goal", {2}}});
}

/*!
 * The states 0 and 1 offer the same two (interval) distributions over the states 1 and 2, so they are bisimilar. State 2 is absorbing and labeled "goal".
 */
std::shared_ptr<storm::models::sparse::Mdp<storm::Interval>> buildIntervalMdp() {
    storm::storage::SparseMatrixBuilder<storm::Interval> builder(5, 3, 0, true, true, 3);
    builder.newRowGroup(0);
    builder.addNextValue(0, 1, storm::Interval(0.1, 0.3));
    builder.addNextValue(0, 2, storm::Interval(0.7, 0.9));
    builder.addNextValue(1, 2, storm::Interval(1.0, 1.0));
    builder.newRowGroup(2);
    builder.addNextValue(2, 1, storm::Interval(0.1, 0.3));
    builder.addNextValue(2, 2, storm::Interval(0.7, 0.9));
    builder.addNextValue(3, 2, storm::Interval(1.0, 1.0));
    builder.newRowGroup(4);
    builder.addNextValue(4, 2, storm::Interval(1.0, 1.0));
    return buildModel<storm::models::sparse::Mdp<storm::Interval>>(builder.build(), {{"goal", {2}}});
}

/*!
 * @return the transition values of the given quotient, in the order of the matrix entries.
 */
template<typename ValueType>
std::vector<ValueType> transitionValues(storm::models::sparse::Model<ValueType> const& model) {
    std::vector<ValueType> result;
    for (auto const& entry : model.getTransitionMatrix()) {
        result.push_back(entry.getValue());
    }
    return result;
}

using Interval = storm::Interval;

/*!
 * Bisimulation on a model whose values already are intervals keeps those intervals.
 */
TEST(IntervalBisimulationTest, IntervalDtmc) {
    auto const model = buildDtmc<Interval>(Interval(0.4, 0.6));
    auto const quotient = storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, strongOptions()).quotient;
    EXPECT_EQ(2ull, quotient->getNumberOfStates());  // {0,1}, {2}
    EXPECT_EQ(std::vector<Interval>({Interval(0.4, 0.6), Interval(0.4, 0.6), Interval(1.0, 1.0)}), transitionValues(*quotient));
}

TEST(IntervalBisimulationTest, IntervalMdp) {
    auto const model = buildIntervalMdp();
    auto const quotient = storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, strongOptions()).quotient;
    EXPECT_EQ(2ull, quotient->getNumberOfStates());  // {0,1}, {2}
    EXPECT_EQ(3ull, quotient->getNumberOfChoices());
}

/*!
 * With a quotient value type that differs from the one of the model, the values of the model are abstracted into intervals.
 * @note the abstraction requires a positive tolerance, as it only pays off if it groups states with (slightly) different values.
 */
TEST(IntervalBisimulationTest, IntervalAbstraction) {
    auto const model = buildDtmc<double>(0.5);
    Options options = strongOptions();
    options.tolerance = storm::utility::convertNumber<storm::RationalNumber>(1e-6);
    auto const quotient = storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, options).quotient;
    EXPECT_EQ(2ull, quotient->getNumberOfStates());  // {0,1}, {2}
    // TODO: As long as the quotient values are taken from the representative state, they only are point intervals, cf. Signatures::extendQuotientData.
    EXPECT_EQ(std::vector<Interval>({Interval(0.5, 0.5), Interval(0.5, 0.5), Interval(1.0, 1.0)}), transitionValues(*quotient));
}

/*!
 * Interval abstraction only pays off with a positive tolerance, as the quotient values are point intervals otherwise.
 */
TEST(IntervalBisimulationTest, IntervalAbstractionWithoutTolerance) {
    auto const model = buildDtmc<double>(0.5);
    STORM_SILENT_EXPECT_THROW((storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, strongOptions())),
                              storm::exceptions::NotSupportedException);
}

/*!
 * Interval values are only supported by signature-based refinement on DTMCs and MDPs.
 */
TEST(IntervalBisimulationTest, UnsupportedBisimulationType) {
    auto const model = buildDtmc<Interval>(Interval(0.4, 0.6));
    STORM_SILENT_EXPECT_THROW(storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, weakOptions()),
                              storm::exceptions::NotSupportedException);
}

}  // namespace

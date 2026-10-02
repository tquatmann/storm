#include "BisimulationTestHelper.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/exceptions/NotSupportedException.h"

namespace {

using storm::test::bisimulation::buildModel;
using storm::test::bisimulation::Options;
using storm::test::bisimulation::strongOptions;
using storm::test::bisimulation::weakOptions;

using Interval = storm::Interval;

/*!
 * Builds a DTMC with one state per entry of distributions: that state moves to the absorbing state labeled "goal" with the first value of its entry and to the
 * absorbing state labeled "sink" with the second one. The "goal" and "sink" states are the last two states of the model.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Dtmc<ValueType>> buildDtmc(std::vector<std::pair<ValueType, ValueType>> const& distributions) {
    uint64_t const goalState = distributions.size();
    uint64_t const sinkState = goalState + 1;
    storm::storage::SparseMatrixBuilder<ValueType> builder(sinkState + 1, sinkState + 1);
    for (uint64_t state = 0; state < distributions.size(); ++state) {
        builder.addNextValue(state, goalState, distributions[state].first);
        builder.addNextValue(state, sinkState, distributions[state].second);
    }
    builder.addNextValue(goalState, goalState, storm::utility::one<ValueType>());
    builder.addNextValue(sinkState, sinkState, storm::utility::one<ValueType>());
    return buildModel<storm::models::sparse::Dtmc<ValueType>>(builder.build(), {{"goal", {goalState}}, {"sink", {sinkState}}});
}

/*!
 * Builds an MDP with one state per entry of distributions: that state offers one choice per distribution of its entry, moving to the absorbing state labeled
 * "goal" with the first value of the distribution and to the absorbing state labeled "sink" with the second one. The "goal" and "sink" states are the last two
 * states of the model.
 */
template<typename ValueType>
std::shared_ptr<storm::models::sparse::Mdp<ValueType>> buildMdp(std::vector<std::vector<std::pair<ValueType, ValueType>>> const& distributions) {
    uint64_t const goalState = distributions.size();
    uint64_t const sinkState = goalState + 1;
    uint64_t numRows = 2;  // the choices of the "goal" and "sink" states
    for (auto const& stateDistributions : distributions) {
        numRows += stateDistributions.size();
    }
    storm::storage::SparseMatrixBuilder<ValueType> builder(numRows, sinkState + 1, 0, true, true, sinkState + 1);
    uint64_t row = 0;
    for (auto const& stateDistributions : distributions) {
        builder.newRowGroup(row);
        for (auto const& [goalValue, sinkValue] : stateDistributions) {
            builder.addNextValue(row, goalState, goalValue);
            builder.addNextValue(row, sinkState, sinkValue);
            ++row;
        }
    }
    builder.newRowGroup(row);
    builder.addNextValue(row++, goalState, storm::utility::one<ValueType>());
    builder.newRowGroup(row);
    builder.addNextValue(row, sinkState, storm::utility::one<ValueType>());
    return buildModel<storm::models::sparse::Mdp<ValueType>>(builder.build(), {{"goal", {goalState}}, {"sink", {sinkState}}});
}

/*!
 * @return the transition values of the given model, in the order of the matrix entries.
 */
template<typename ValueType>
std::vector<ValueType> transitionValues(storm::models::sparse::Model<ValueType> const& model) {
    std::vector<ValueType> result;
    for (auto const& entry : model.getTransitionMatrix()) {
        result.push_back(entry.getValue());
    }
    return result;
}

/*!
 * Checks that the transition values of the given model are the given intervals. The bounds are compared up to a small numerical tolerance, since restricting
 * an interval to the feasible values can move a bound by a few ulps.
 */
void expectTransitionValues(storm::models::sparse::Model<Interval> const& model, std::vector<Interval> const& expected) {
    auto const actual = transitionValues(model);
    ASSERT_EQ(expected.size(), actual.size());
    for (uint64_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i].lower(), actual[i].lower(), 1e-12) << "lower bound of transition " << i;
        EXPECT_NEAR(expected[i].upper(), actual[i].upper(), 1e-12) << "upper bound of transition " << i;
    }
}

/*!
 * @return options with the given (positive) tolerance, which is what the abstraction of values into intervals is allowed to deviate by.
 */
Options optionsWithTolerance(double const tolerance) {
    Options options = strongOptions();
    options.tolerance = storm::utility::convertNumber<storm::RationalNumber>(tolerance);
    return options;
}

storm::RationalNumber rational(std::string const& number) {
    return storm::utility::convertNumber<storm::RationalNumber>(number);
}

// ------------------------------------------------------------
// Models that already have interval values
// ------------------------------------------------------------

/*!
 * Bisimulation on a model whose values already are intervals keeps those intervals.
 */
TEST(IntervalBisimulationTest, IntervalDtmc) {
    auto const model = buildDtmc<Interval>({{Interval(0.4, 0.6), Interval(0.4, 0.6)}, {Interval(0.4, 0.6), Interval(0.4, 0.6)}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, strongOptions()).quotient;
    EXPECT_EQ(3ull, quotient->getNumberOfStates());  // {0,1}, {goal}, {sink}
    expectTransitionValues(*quotient, {Interval(0.4, 0.6), Interval(0.4, 0.6), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

/*!
 * The two states offer the same two interval distributions, so they are bisimilar. Their choices are kept apart, as the quotient may not offer less behaviour.
 */
TEST(IntervalBisimulationTest, IntervalMdp) {
    auto const model = buildMdp<Interval>({{{Interval(0.1, 0.3), Interval(0.7, 0.9)}, {Interval(1.0, 1.0), Interval(0.0, 0.0)}},
                                           {{Interval(0.1, 0.3), Interval(0.7, 0.9)}, {Interval(1.0, 1.0), Interval(0.0, 0.0)}}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, strongOptions()).quotient;
    EXPECT_EQ(3ull, quotient->getNumberOfStates());  // {0,1}, {goal}, {sink}
    EXPECT_EQ(4ull, quotient->getNumberOfChoices());
}

/*!
 * The bounds of an interval distribution are restricted to the values that can actually be instantiated to a probability distribution: the first state can move
 * to "goal" with at most 1 - 0.5 and to "sink" with at most 1 - 0.4, which makes it behave exactly like the second state.
 */
TEST(IntervalBisimulationTest, InfeasibleIntervalsAreRestricted) {
    auto const model = buildDtmc<Interval>({{Interval(0.4, 0.9), Interval(0.5, 0.8)}, {Interval(0.4, 0.5), Interval(0.5, 0.6)}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, strongOptions()).quotient;
    EXPECT_EQ(3ull, quotient->getNumberOfStates());  // {0,1}, {goal}, {sink}
    expectTransitionValues(*quotient, {Interval(0.4, 0.5), Interval(0.5, 0.6), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

// ------------------------------------------------------------
// Abstraction of the values of a model into intervals
// ------------------------------------------------------------

/*!
 * With a quotient value type that differs from the one of the model, the values of the model are abstracted into intervals: the two states are grouped
 * although their values differ, and the quotient covers the values of both of them.
 */
TEST(IntervalBisimulationTest, AbstractionWidensToConvexHull) {
    auto const model = buildDtmc<double>({{0.5, 0.5}, {0.5008, 0.4992}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, optionsWithTolerance(0.002)).quotient;
    EXPECT_EQ(3ull, quotient->getNumberOfStates());  // {0,1}, {goal}, {sink}
    expectTransitionValues(*quotient, {Interval(0.5, 0.5008), Interval(0.4992, 0.5), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

/*!
 * Grouping states transitively would widen the intervals without a bound, so the amount a signature is widened by is limited: the first two states are grouped
 * (0.0008 apart), but the third one is not, as covering it as well would widen the interval by more than the tolerance allows.
 */
TEST(IntervalBisimulationTest, AbstractionIsBoundedByTheTolerance) {
    auto const model = buildDtmc<double>({{0.5, 0.5}, {0.5008, 0.4992}, {0.5016, 0.4984}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, optionsWithTolerance(0.002)).quotient;
    EXPECT_EQ(4ull, quotient->getNumberOfStates());  // {0,1}, {2}, {goal}, {sink}
    expectTransitionValues(
        *quotient, {Interval(0.5, 0.5008), Interval(0.4992, 0.5), Interval(0.5016, 0.5016), Interval(0.4984, 0.4984), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

/*!
 * States whose values are further apart than the tolerance are not grouped, so their values are not abstracted at all.
 */
TEST(IntervalBisimulationTest, AbstractionKeepsDistantStatesApart) {
    auto const model = buildDtmc<double>({{0.5, 0.5}, {0.505, 0.495}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, optionsWithTolerance(0.002)).quotient;
    EXPECT_EQ(4ull, quotient->getNumberOfStates());  // {0}, {1}, {goal}, {sink}
    expectTransitionValues(*quotient,
                           {Interval(0.5, 0.5), Interval(0.5, 0.5), Interval(0.505, 0.505), Interval(0.495, 0.495), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

/*!
 * Regression test: restricting the point intervals of this state to the feasible values yields an empty intersection, as 1 - 0.4984 is a few ulps below 0.5016.
 * An interval whose lower bound exceeds its upper bound has NaN bounds, so that case has to be caught before the bounds are written.
 */
TEST(IntervalBisimulationTest, AbstractionOfInfeasibleRoundedValues) {
    auto const model = buildDtmc<double>({{0.5016, 0.4984}});
    auto const quotient = storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, optionsWithTolerance(1e-6)).quotient;
    EXPECT_EQ(3ull, quotient->getNumberOfStates());  // nothing to group here
    expectTransitionValues(*quotient, {Interval(0.5016, 0.5016), Interval(0.4984, 0.4984), Interval(1.0, 1.0), Interval(1.0, 1.0)});
}

/*!
 * The values of an exact model are abstracted into rational intervals.
 */
TEST(IntervalBisimulationTest, AbstractionOfExactValues) {
    auto const model = buildDtmc<storm::RationalNumber>(
        {{rational("1/2"), rational("1/2")}, {rational("501/1000"), rational("499/1000")}, {rational("3/4"), rational("1/4")}});
    auto const quotient =
        storm::bisimulation::performBisimulationMinimization<storm::RationalNumber, storm::RationalInterval>(*model, {}, optionsWithTolerance(0.01)).quotient;
    EXPECT_EQ(4ull, quotient->getNumberOfStates());  // {0,1}, {2}, {goal}, {sink}
    std::vector<storm::RationalInterval> const expected{
        storm::RationalInterval(rational("1/2"), rational("501/1000")), storm::RationalInterval(rational("499/1000"), rational("1/2")),
        storm::RationalInterval(rational("3/4"), rational("3/4")),      storm::RationalInterval(rational("1/4"), rational("1/4")),
        storm::RationalInterval(rational("1"), rational("1")),          storm::RationalInterval(rational("1"), rational("1"))};
    EXPECT_EQ(expected, transitionValues(*quotient));
}

// ------------------------------------------------------------
// Unsupported configurations
// ------------------------------------------------------------

/*!
 * Abstracting values into intervals only pays off with a positive tolerance, as the quotient values are point intervals otherwise.
 */
TEST(IntervalBisimulationTest, AbstractionRequiresTolerance) {
    auto const model = buildDtmc<double>({{0.5, 0.5}, {0.5, 0.5}});
    STORM_SILENT_EXPECT_THROW((storm::bisimulation::performBisimulationMinimization<double, Interval>(*model, {}, strongOptions())),
                              storm::exceptions::NotSupportedException);
}

/*!
 * Interval values are only supported by signature-based refinement on DTMCs and MDPs.
 */
TEST(IntervalBisimulationTest, UnsupportedBisimulationType) {
    auto const model = buildDtmc<Interval>({{Interval(0.4, 0.6), Interval(0.4, 0.6)}, {Interval(0.4, 0.6), Interval(0.4, 0.6)}});
    STORM_SILENT_EXPECT_THROW(storm::bisimulation::performBisimulationMinimization<Interval>(*model, {}, weakOptions()),
                              storm::exceptions::NotSupportedException);
}

}  // namespace

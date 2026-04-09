#include <sstream>
#include "storm/models/sparse/Dtmc.h"
#include "storm/storage/bisimulation/Partition.h"

#include "test/storm_gtest.h"

#include <storm-parsers/parser/DirectEncodingParser.h>
#include <memory>
#include <vector>
#include "storm/adapters/IntervalAdapter.h"
#include "storm/storage/geometry/Halfspace.h"
#include "storm/storage/geometry/Polytope.h"

namespace {

std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>> create2DPolytope(storm::RationalNumber c1LowerBound,
                                                                                            storm::RationalNumber c1UpperBound,
                                                                                            storm::RationalNumber c2LowerBound,
                                                                                            storm::RationalNumber c2UpperBound) {
    using Point = typename storm::storage::geometry::Polytope<storm::RationalNumber>::Point;

    // Halfspace in R² is given by: a1 * p1 + a2 * p2 <= b, where (a1, a2) is the normal vector and b is the offset
    std::vector<storm::storage::geometry::Halfspace<storm::RationalNumber>> halfspaces;

    // Interval bounds
    halfspaces.emplace_back(Point{-1.0, 0.0},  // a1, a2
                            -c1LowerBound);    // => (-p1 <= -l1) <=> (p1 >= l1)
    halfspaces.emplace_back(Point{1.0, 0.0},   // a1, a2
                            c1UpperBound);     // => (p1 <= u1)
    halfspaces.emplace_back(Point{0.0, -1.0},  // a1, a2
                            -c2LowerBound);    // => (-p2 <= -l2) <=> (p2 >= l2)
    halfspaces.emplace_back(Point{0.0, 1.0},   // a1, a2
                            c2UpperBound);     // => (p2 <= u2)

    // Normalization constraint: p1 + p2 = 1
    halfspaces.emplace_back(Point{1.0, 1.0},    // a1, a2
                            1.0);               // p1 + p2 <= 1
    halfspaces.emplace_back(Point{-1.0, -1.0},  // a1, a2
                            -1.0);              // -(p1 + p2) <= -1  => p1 + p2 >= 1

    // Polytope is "intersection" of all halfspaces
    return storm::storage::geometry::Polytope<storm::RationalNumber>::create(halfspaces);
}

void print2DPolytope(const std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>>& polytope) {
    ASSERT_FALSE(polytope->isEmpty());
    ASSERT_EQ(2, polytope->getVertices().size());

    auto i = 1;
    for (auto vertice : polytope->getVertices()) {
        std::cout << "Vertice " << i << ": " << "(" << vertice[0].get_d() << ", " << vertice[1].get_d() << ")" << std::endl;
        i++;
    }
}

TEST(PolytopeTest, IDTMC2DPolytope) {
    // Consider the following two states s_1, s_2, where
    // I(s_1, C_1) = [0.1, 0.5], I(s_1, C_2) = [0.8, 1]
    // I(s_2, C_1) = [0.1, 0.3], I(s_1, C_2) = [0.7, 1]

    auto polytopeState1 = create2DPolytope(0.1, 0.5, 0.8, 1.0);
    auto polytopeState2 = create2DPolytope(0.1, 0.3, 0.7, 1.0);

    ASSERT_FALSE(polytopeState1->isEmpty());
    ASSERT_FALSE(polytopeState2->isEmpty());

    std::cout << "Polytope of state 1: " << std::endl;
    print2DPolytope(polytopeState1);

    std::cout << "Polytope of state 2: " << std::endl;
    print2DPolytope(polytopeState2);

    ASSERT_FALSE(polytopeState1->contains(polytopeState2));
    ASSERT_TRUE(polytopeState2->contains(polytopeState1));

    auto intersectionOfStatePolytopes = polytopeState1->intersection(polytopeState2);
    std::cout << "Polytope of intersection: " << std::endl;
    print2DPolytope(intersectionOfStatePolytopes);

    ASSERT_TRUE(intersectionOfStatePolytopes->contains(polytopeState1));
    ASSERT_TRUE(polytopeState1->contains(intersectionOfStatePolytopes));
}

TEST(PolytopeTest, CreatePolytopesfromIDTMC) {
    // Just tests the parsing of an IDTMC and creates useless polytopes from it.

    std::shared_ptr<storm::models::sparse::Model<storm::Interval>> modelPtr =
        storm::parser::parseDirectEncodingModel<storm::Interval>(STORM_TEST_RESOURCES_DIR "/idtmc/brp-16-2.drn");
    std::shared_ptr<storm::models::sparse::Dtmc<storm::Interval>> dtmc = modelPtr->as<storm::models::sparse::Dtmc<storm::Interval>>();
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(613ul, dtmc->getNumberOfStates());
    EXPECT_TRUE(modelPtr->hasUncertainty());

    std::vector<std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>>> polytopesOfStates(dtmc->getNumberOfStates());

    for (auto i = 0; i < dtmc->getNumberOfStates(); i++) {
        for (const auto& entry : dtmc->getTransitionMatrix().getRow(i)) {
            polytopesOfStates.emplace_back(create2DPolytope(entry.getValue().lower(), entry.getValue().upper(), 0.0, 1.0));
        }
    }

    std::cout << polytopesOfStates.size() << std::endl;
}

}  // namespace
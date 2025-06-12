#include "DeterministicIntervalModelBisimulationDecomposition.h"
#include <algorithm>
#include <boost/iterator/zip_iterator.hpp>
#include <chrono>
#include <iomanip>
#include <storm/models/sparse/Dtmc.h>
#include <unordered_map>
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"

#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/graph.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::DeterministicIntervalModelBisimulationDecomposition(
    const ModelType& model, const typename BisimulationDecomposition<ModelType>::Options& options)
    : storm::storage::BisimulationDecomposition<ModelType>(model, options) {}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::buildQuotient() {
    std::cout << "Size of partition: " << this->partition.getNumberOfBlocks();
    // TODO: implement actual quotient model construction
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::initializeMeasureDrivenPartition() {
    // TODO: implement actual measure-driven partitioning
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::initializeLabelBasedPartition() {
    BisimulationDecomposition<ModelType>::initializeLabelBasedPartition();
    postProcessInitialPartition();
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::postProcessInitialPartition() {
    // TODO: implement
}

template<typename ModelType>
void DeterministicIntervalModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(
    std::span<uint64_t const> splitterBlock,
    std::deque<typename bisimulation::Partition::Block>& splitterQueue,
    bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::bisimulation::Partition::BlockSet blocksToSplit;

    std::vector<double> lowerBoundOfStatesToSplitter(this->model.getNumberOfStates());
    std::vector<double> upperBoundOfStatesToSplitter(this->model.getNumberOfStates());

    std::vector<double> lowerBoundOfStatesToOtherBlocks(this->model.getNumberOfStates());
    std::vector<double> upperBoundOfStatesToOtherBlocks(this->model.getNumberOfStates());

    std::vector<std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>>> polytopesOfStates(this->model.getNumberOfStates());

    // get all predecessor blocks
    for (auto currentState : splitterBlock) {
        // Compute probability to enter splitter block for each predecessor
        for (const auto& predecessorEntry : this->backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);

            lowerBoundOfStatesToSplitter[predecessorEntry.getColumn()] += predecessorEntry.getValue().lower();
            upperBoundOfStatesToSplitter[predecessorEntry.getColumn()] += predecessorEntry.getValue().upper();

            // Remember which blocks contain predecessors to split them w.r.t. the splitter afterwards
            blocksToSplit.insert(predecessorBlock);
        }
    }

    // Loop through all predecessor blocks
    for (auto predecessorBlockToSplit : blocksToSplit) {
        // Loop through all predecessor states
        for (auto currentState : predecessorBlockToSplit) {
            // Loop through all successors of a predecessor
            for (const auto& entry : this->model.getTransitionMatrix().getRow(currentState)) {
                auto successorState = entry.getColumn();
                auto blockOfSuccessor = this->partition.getBlockOfElement(successorState);

                // Transition to another block and NOT the splitter block
                if (!this->partition.contains(splitterBlock, successorState)) {
                    lowerBoundOfStatesToOtherBlocks[currentState] += entry.getValue().lower();
                    upperBoundOfStatesToOtherBlocks[currentState] += entry.getValue().upper();
                }
            }
        }
    }

    // Loop through all predecessor blocks
    for (auto predecessorBlockToSplit : blocksToSplit) {
        // Loop through all predecessor states
        for (auto currentState : predecessorBlockToSplit) {
            polytopesOfStates[currentState] = create2DPolytope(std::min(1.0, lowerBoundOfStatesToSplitter[currentState]), std::min(1.0, upperBoundOfStatesToSplitter[currentState]),
                             std::min(1.0, lowerBoundOfStatesToOtherBlocks[currentState]), std::min(1.0, upperBoundOfStatesToOtherBlocks[currentState]));
        }
    }

    // Loop through all predecessor blocks
    for (auto predecessorBlockToSplit : blocksToSplit) {
        std::unordered_set<uint64_t> equivalentPolytopes;
        auto stateToCompare = predecessorBlockToSplit.front();
        equivalentPolytopes.insert(stateToCompare);
        for (auto predecessorState : predecessorBlockToSplit) {
            if (predecessorState == stateToCompare) {
                continue;
            }

            if (polytopesOfStates[stateToCompare]->contains(polytopesOfStates[predecessorState]) &&
                polytopesOfStates[predecessorState]->contains(polytopesOfStates[stateToCompare])) {
                equivalentPolytopes.insert(predecessorState);
            }
        }

        auto [notEquivalentBlock, equivalentBlock] =
            this->partition.splitBlockByPredicate(predecessorBlockToSplit, [&equivalentPolytopes]
                                                  (auto const& state) { return equivalentPolytopes.contains(state); });

        if (notEquivalentBlock.size() > 0 && !enqueuedSplitterBlocks.contains(notEquivalentBlock)) {
            enqueuedSplitterBlocks.insert(notEquivalentBlock);
            splitterQueue.push_back(notEquivalentBlock);
        }

        if (equivalentBlock.size() > 0 && !enqueuedSplitterBlocks.contains(equivalentBlock)) {
            enqueuedSplitterBlocks.insert(equivalentBlock);
            splitterQueue.push_back(equivalentBlock);
        }
    }
}

template<typename ModelType>
std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>>
    DeterministicIntervalModelBisimulationDecomposition<ModelType>::create2DPolytope(storm::RationalNumber c1LowerBound, storm::RationalNumber c1UpperBound,
                                                                                            storm::RationalNumber c2LowerBound, storm::RationalNumber c2UpperBound) {
    using Point = typename storm::storage::geometry::Polytope<storm::RationalNumber>::Point;

    // Halfspace in R² is given by: a1 * p1 + a2 * p2 <= b, where (a1, a2) is the normal vector and b is the offset
    std::vector<storm::storage::geometry::Halfspace<storm::RationalNumber>> halfspaces;

    // Interval bounds
    halfspaces.emplace_back(Point{-1.0, 0.0}, // a1, a2
                            -c1LowerBound);  // => (-p1 <= -l1) <=> (p1 >= l1)
    halfspaces.emplace_back(Point{ 1.0, 0.0}, // a1, a2
                            c1UpperBound);  // => (p1 <= u1)
    halfspaces.emplace_back(Point{ 0.0,-1.0}, // a1, a2
                            -c2LowerBound);  // => (-p2 <= -l2) <=> (p2 >= l2)
    halfspaces.emplace_back(Point{ 0.0, 1.0}, // a1, a2
                            c2UpperBound);  // => (p2 <= u2)

    // Normalization constraint: p1 + p2 = 1
    halfspaces.emplace_back(Point{ 1.0, 1.0}, // a1, a2
                            1.0);  // p1 + p2 <= 1
    halfspaces.emplace_back(Point{-1.0,-1.0}, // a1, a2
                            -1.0);  // -(p1 + p2) <= -1  => p1 + p2 >= 1

    // Polytope is "intersection" of all halfspaces
    return storm::storage::geometry::Polytope<storm::RationalNumber>::create(halfspaces);
}

template<typename ModelType>
std::pair<storm::storage::BitVector, storm::storage::BitVector>
DeterministicIntervalModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    // TODO: return actual state sets for probability 0 and 1
    return {storm::storage::BitVector(), storm::storage::BitVector()};
}

template class DeterministicIntervalModelBisimulationDecomposition<storm::models::sparse::Dtmc<storm::Interval>>;

}  // namespace storage
}  // namespace storm

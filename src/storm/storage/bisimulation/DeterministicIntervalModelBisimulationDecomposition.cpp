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
    : storm::storage::BisimulationDecomposition<ModelType>(model, options),
      probabilitiesToCurrentSplitter(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      probabilitiesToOtherBlocks(model.getNumberOfStates(), storm::utility::zero<ValueType>()),
      touchedProbabilitiesToSplitter(model.getNumberOfStates(), false) {}

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
    storm::storage::Partition::Block splitterBlock,
    std::deque<typename bisimulation::Partition::Block>& splitterQueue,
    bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) {
    storm::storage::bisimulation::Partition::BlockSet blocksToSplit;

    // std::fill(probabilitiesToCurrentSplitter.begin(), probabilitiesToCurrentSplitter.end(), storm::utility::zero<ValueType>());
    for (auto currentState : splitterBlock) {
        // Compute probability to enter splitter block for each predecessor
        for (const auto& predecessorEntry : this->backwardTransitions.getRow(currentState)) {
            auto predecessorState = predecessorEntry.getColumn();
            auto predecessorBlock = this->partition.getBlockOfElement(predecessorState);
            auto intervalTransitionProbability = predecessorEntry.getValue();

            if (!possiblyNeedsRefinement(predecessorBlock)) {
                continue;
            }

            // TODO: Compute both intervals going to the splitter and going to all other blocks
            if (touchedProbabilitiesToSplitter.get(predecessorState)) {
                probabilitiesToCurrentSplitter[predecessorState] += intervalTransitionProbability;
            } else {
                probabilitiesToCurrentSplitter[predecessorState] = intervalTransitionProbability;
                probabilitiesToOtherBlocks[predecessorState] = storm::utility::zero<ValueType>(); // init
                touchedProbabilitiesToSplitter.set(predecessorState, true);
            }

            // Remember which blocks contain predecessors to split them w.r.t. the splitter afterwards
            blocksToSplit.insert(predecessorBlock);
        }
    }

    for (auto predecessorState : touchedProbabilitiesToSplitter) {
        for (const auto& entry : this->model.getTransitionMatrix().getRow(predecessorState)) {
            auto targetState = entry.getColumn();

            if (!this->partition.contains(splitterBlock, targetState)) {
                probabilitiesToOtherBlocks[predecessorState] += entry.getValue(); // interval addition
            }
        }
    }

    //for (auto predecessorBlockToSplit : blocksToSplit) {
    //    // First split the block by whether it is a predecessor of the splitter block or not
    //    bool wasSplit = this->partition.splitBlockByOrder(predecessorBlockToSplit, [this]
    //                                                      (auto const& a, auto const& b) {
    //                                                          if (!touchedProbabilitiesToSplitter.get(a) || !touchedProbabilitiesToSplitter.get(b)) {
    //                                                              std::cout << "Comparing states that were not touched: " << a << ", " << b << std::endl;
    //                                                              return false;
    //                                                          }
    //                                                          return computeIntervalProjection(probabilitiesToCurrentSplitter[a], probabilitiesToOtherBlocks[a])
    //                                                          < computeIntervalProjection(probabilitiesToCurrentSplitter[b], probabilitiesToOtherBlocks[b]);
    //                                                      });
//
    //    // Add all blocks that were split to splitter queue
    //    if (wasSplit) {
    //        this->partition.forEachSubBlock(predecessorBlockToSplit, [&splitterQueue, &enqueuedSplitterBlocks](auto const& block) {
    //            if (!enqueuedSplitterBlocks.contains(block)) {
    //                splitterQueue.push_back(block);
    //                enqueuedSplitterBlocks.insert(block);
    //            }
    //        });
    //    }
    //}

    for (auto predecessorBlockToSplit : blocksToSplit) {
        // First split the block by whether it is a predecessor of the splitter block or not
        auto [noPredecessors, predecessors] =
            this->partition.splitBlockByPredicate(predecessorBlockToSplit, [this]
                                                  (auto const& state) { return touchedProbabilitiesToSplitter.get(state); });

        if (noPredecessors.size() > 0) {
            if (!enqueuedSplitterBlocks.contains(noPredecessors)) {
                splitterQueue.push_back(noPredecessors);
                enqueuedSplitterBlocks.insert(noPredecessors);
            }
        }

        if (predecessors.size() > 0) {
            bool wasSplit = this->partition.splitBlockByOrder(predecessors, [this]
                                                              (auto const& a, auto const& b) {
                                                                  if (!touchedProbabilitiesToSplitter.get(a) || !touchedProbabilitiesToSplitter.get(b)) {
                                                                      std::cout << "Comparing states that were not touched: " << a << ", " << b << std::endl;
                                                                      return false;
                                                                  }

                                                                  return computeIntervalProjection(probabilitiesToCurrentSplitter[a], probabilitiesToOtherBlocks[a])
                                                                         < computeIntervalProjection(probabilitiesToCurrentSplitter[b], probabilitiesToOtherBlocks[b]);
                                                              });

            // Add all blocks that were split to splitter queue
            if (wasSplit) {
                this->partition.forEachSubBlock(predecessors, [&splitterQueue, &enqueuedSplitterBlocks](auto const& block) {
                    if (!enqueuedSplitterBlocks.contains(block)) {
                        splitterQueue.push_back(block);
                        enqueuedSplitterBlocks.insert(block);
                    }
                });
            }
        }
    }

    // Reset the touched entries of the probabilitiesToCurrentSplitter vector
    touchedProbabilitiesToSplitter.clear();
}

template<typename ModelType>
ModelType::ValueType DeterministicIntervalModelBisimulationDecomposition<ModelType>::computeIntervalProjection(ValueType intervalToSplitter, ValueType intervalToOtherBlocks) {
    // Normalize intervals
    storm::Interval normalizedIntervalToSplitter(std::min(1.0, intervalToSplitter.lower()), std::min(1.0, intervalToSplitter.upper()));
    storm::Interval normalizedIntervalToOtherBlocks(std::min(1.0, intervalToOtherBlocks.lower()), std::min(1.0, intervalToOtherBlocks.upper()));

    // Compute interval projection
    return storm::Interval(std::max(normalizedIntervalToSplitter.lower(), 1.0 - normalizedIntervalToOtherBlocks.upper()),
                           std::min(normalizedIntervalToSplitter.upper(), 1.0 - normalizedIntervalToOtherBlocks.lower()));
}

template<typename ModelType>
bool DeterministicIntervalModelBisimulationDecomposition<ModelType>::possiblyNeedsRefinement(std::span<uint64_t const> block) const {
    return block.size() > 1 && !this->absorbingBlocks.contains(block.front());
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

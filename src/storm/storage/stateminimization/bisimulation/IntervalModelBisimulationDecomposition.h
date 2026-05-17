#pragma once

#include <boost/iterator/zip_iterator.hpp>
#include <chrono>
#include <unordered_map>

#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionForward.h"
#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/DeterministicModel.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/Distribution.h"
#include "storm/storage/geometry/Polytope.h"
#include "storm/storage/stateminimization/Partition.h"
#include "storm/storage/stateminimization/bisimulation/BisimulationDecomposition.h"
#include "storm/utility/constants.h"
#include "storm/utility/graph.h"

namespace storm {
namespace storage {

// TODO: Extend logic to support deterministic and nondeterministic models.
template<typename ModelType>
class IntervalModelBisimulationDecomposition : public BisimulationDecomposition<ModelType> {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef typename ModelType::RewardModelType RewardModelType;

    /*!
     * Computes the bisimulation relation for the given model. Which kind of bisimulation is computed, is
     * customizable via the given options.
     *
     * @param model The model to decompose.
     * @param options The options that customize the computed bisimulation.
     */
    IntervalModelBisimulationDecomposition(ModelType const& model, typename BisimulationDecomposition<ModelType>::BisimulationOptions const& options =
                                                                       typename BisimulationDecomposition<ModelType>::BisimulationOptions());

   protected:
    std::pair<storm::storage::BitVector, storm::storage::BitVector> getStatesWithProbability01() override;

    void refinePartitionBasedOnSplitter(std::span<uint64_t const> splitterBlock, std::deque<typename stateminimization::Partition::Block>& splitterQueue,
                                        stateminimization::Partition::BlockSet& enqueuedSplitterBlocks) override;

    void buildQuotientFromPartition() override;

   private:
    // Retrieves whether the given predecessor of the splitters possibly needs refinement.
    bool possiblyNeedsRefinement(std::span<uint64_t const> block) const;

    /*!
     * Computes the feasible interval based on the intervals going to the splitter block C and all other blocks \Pi \setminus C.
     * @param intervalToSplitter aggregated interval to C
     * @param intervalToOtherBlocks aggregated interval \Pi \setminus C
     * @return feasible interval to splitter block C
     */
    ValueType computeFeasibleIntervalBasedOnAggregatedIntervals(ValueType intervalToSplitter, ValueType intervalToOtherBlocks) const;

    /*!
     * Takes the given interval and computes the intersection with [0, 1] with respect to the given precision.
     * @param interval
     * @return probabilistic interval I s.t. I = interval \cap [0, 1]
     */
    ValueType clampIntervalToProbabilisticInterval(ValueType interval) const;

    /*!
     * Returns the given distribution where each interval entry I gets clamped to a probabilistic interval I', i.e., I' = I \cap [0, 1].
     * @param distribution
     * @return distribution with clamped interval entries
     */
    Distribution<ValueType, sparse::state_type> getClampedDistribution(Distribution<ValueType, storm::storage::sparse::state_type> distribution) const;

    bool checkCurrentPartitionByExactFeasibleIntervals() const;

    /// A vector that holds the probabilities of states going into the splitter. This is used by the method that
    /// refines a block based on probabilities.
    std::vector<ValueType> probabilitiesToCurrentSplitter;

    /// A vector that holds the probabilities of states going into the splitter. This is used by the method that
    /// refines a block based on probabilities.
    std::vector<ValueType> probabilitiesToOtherBlocks;

    /// A bitvector indicating which positions of the probabilitiesToCurrentSplitter were touched in this iteration.
    /// This is an alternative solution to the marker1-logic that was part of the old implementation.
    /// Note that this bitvector also indicates whether a state is a direct predecessor of the splitter block or not.
    storm::storage::BitVector touchedProbabilitiesToSplitter;

    /// A vector that maps each choice entry to its owner state.
    std::vector<storm::storage::sparse::state_type> choiceToStateMapping;

    /// The backward transitions of the model without merging the choice rows into state representation.
    storm::storage::SparseMatrix<ValueType> backwardTransitionsWithChoices;
};

}  // namespace storage
}  // namespace storm

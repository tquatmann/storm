#ifndef STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H
#define STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H

#include <storm/storage/Distribution.h>
#include <storm/storage/bisimulation/Partition.h>
#include <storm/storage/geometry/Polytope.h>
#include "storm/adapters/IntervalForward.h"
#include "storm/models/sparse/DeterministicModel.h"
#include "storm/storage/bisimulation/BisimulationDecomposition.h"

namespace storm {
namespace storage {

template<typename ModelType>
class DeterministicIntervalModelBisimulationDecomposition : public BisimulationDecomposition<ModelType> {
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
    DeterministicIntervalModelBisimulationDecomposition(ModelType const& model, typename BisimulationDecomposition<ModelType>::Options const& options =
                                                                                    typename BisimulationDecomposition<ModelType>::Options());

   protected:
    virtual std::pair<storm::storage::BitVector, storm::storage::BitVector> getStatesWithProbability01() override;

    virtual void initializeMeasureDrivenPartition() override;

    virtual void initializeLabelBasedPartition() override;

    virtual void buildQuotient() override;

    virtual void refinePartitionBasedOnSplitter(std::span<uint64_t const> splitterBlock, std::deque<typename bisimulation::Partition::Block>& splitterQueue,
                                                bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) override;

    virtual void refineBlockBasedOnEpsilonSignature(std::span<uint64_t const> subBlock, std::deque<typename bisimulation::Partition::Block>& blocksQueue,
                                                    bisimulation::Partition::BlockSet& enqueuedBlocks, double epsilon) override;

   private:
    void postProcessInitialPartition();

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

    /*!
     * Computes the delta between the distribution of a state and its candidate group distribution.
     * @param state state of interest
     * @param stateDistribution distribution over the partition of the state
     * @param enhancedDistribution enhanced distribution over the partition based on its (candidate) group
     * @return delta
     */
    storm::IntervalBaseType<ValueType> computeDeltaForState(
        storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& stateDistribution,
        storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& enhancedDistribution);

    /*!
     * Takes the two input distributions and constructs an enhanced version, i.e., for each entry 'block' it chooses:
     * [ min(\inf firstDistribution(block), \inf secondDistribution(block)), max(\sup firstDistribution(block), \sup secondDistribution(block)) ]
     * @param firstDistribution
     * @param secondDistribution
     * @return enhanced distribution
     */
    storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> computeEnhancedDistribution(
        storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& firstDistribution,
        storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const& secondDistribution);

    bool canMergeBlocksForDebug(std::span<uint64_t const> blockA, std::span<uint64_t const> blockB, double epsilon);

    void debugFindMergeableFinalBlocks(double epsilon);

    void validateEpsilonStability(double epsilon);

    bool checkCurrentPartitionByExactFeasibleIntervals() const;

    // A vector that holds the probabilities of states going into the splitter. This is used by the method that
    // refines a block based on probabilities.
    std::vector<ValueType> probabilitiesToCurrentSplitter;

    // A vector that holds the probabilities of states going into the splitter. This is used by the method that
    // refines a block based on probabilities.
    std::vector<ValueType> probabilitiesToOtherBlocks;

    // A bitvector indicating which positions of the probabilitiesToCurrentSplitter were touched in this iteration.
    // This is an alternative solution to the marker1-logic that was part of the old implementation.
    // Note that this bitvector also indicates whether a state is a direct predecessor of the splitter block or not.
    storm::storage::BitVector touchedProbabilitiesToSplitter;

    // A vector of the distributions over the current partition of each state
    std::vector<storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>> originalStateDistributions;
};

}  // namespace storage
}  // namespace storm

#endif  // STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H

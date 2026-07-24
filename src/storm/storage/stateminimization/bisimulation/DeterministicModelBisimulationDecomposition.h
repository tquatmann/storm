#pragma once

#include <boost/iterator/zip_iterator.hpp>
#include <chrono>
#include <deque>
#include <span>
#include <unordered_map>
#include <map>

#include <boost/container/flat_map.hpp>

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/stateminimization/Partition.h"
#include "storm/storage/stateminimization/bisimulation/BisimulationDecomposition.h"
#include "storm/utility/constants.h"
#include "storm/utility/graph.h"

namespace storm {
namespace storage {

/*!
 * This class represents the decomposition of a deterministic model into its bisimulation quotient.
 */
template<typename ModelType>
class DeterministicModelBisimulationDecomposition : public BisimulationDecomposition<ModelType> {
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
    DeterministicModelBisimulationDecomposition(ModelType const& model, typename BisimulationDecomposition<ModelType>::BisimulationOptions const& options =
                                                                            typename BisimulationDecomposition<ModelType>::BisimulationOptions());

   protected:
    virtual std::pair<storm::storage::BitVector, storm::storage::BitVector> getStatesWithProbability01() override;

    virtual void buildQuotientFromPartition() override;

    virtual void refinePartitionBasedOnSplitter(std::span<uint64_t const> splitterBlock,
                                                std::deque<typename stateminimization::Partition::Block>& splitterQueue,
                                                stateminimization::Partition::OrderedBlockSet& enqueuedSplitterBlocks) override;

   private:

    struct RefinementCache {
        std::vector<ValueType> probabilitiesToSplitter;
        std::vector<uint64_t> splitterPredecessors; // states with a non-zero probability to a splitter
        storm::storage::stateminimization::Partition::NonSuperBlockSet nonSuperBlockSet;

        RefinementCache(storm::storage::stateminimization::Partition const& partition) : probabilitiesToSplitter(partition.getNumberOfElements(), storm::utility::zero<ValueType>()), nonSuperBlockSet(partition) {}

        void addProbabilityToSplitter(uint64_t state, ValueType const& probability) {
            STORM_LOG_ASSERT(!storm::utility::isZero(probability), "The probability to add to the splitter must not be zero.");
            if (storm::utility::isZero(probabilitiesToSplitter[state])) {
                splitterPredecessors.push_back(state);
            }
            probabilitiesToSplitter[state] += probability;
        }

        void clear() {
            for (auto const& state : splitterPredecessors) {
                probabilitiesToSplitter[state] = storm::utility::zero<ValueType>();
            }
            splitterPredecessors.clear();
            STORM_LOG_ASSERT(std::all_of(probabilitiesToSplitter.begin(), probabilitiesToSplitter.end(), [](ValueType const& p) { return storm::utility::isZero(p); }), "Expected all probabilities to be zero after clearing the cache.");
            while (!nonSuperBlockSet.empty()) {
                nonSuperBlockSet.pop();
            }
        }


    } refinementCache;

    /*!
     * Performs the necessary steps to compute a weak bisimulation on a DTMC.
     */
    void initializeWeakDtmcBisimulation();

    /*!
     * Splits all blocks of the current partition into a block that contains all divergent states and another
     * block containing the non-divergent states.
     */
    void splitOffDivergentStates();

    /*!
     * Initializes the vector of silent probabilities.
     */
    void initializeSilentProbabilities();

    // Retrieves whether the given predecessor of the splitters possibly needs refinement.
    bool possiblyNeedsRefinement(std::span<uint64_t const> block) const;

    ValueType getTransitionValue(storm::storage::MatrixEntry<storm::storage::sparse::state_type, ValueType> const& matrixEntry,
                                 [[maybe_unused]] storm::storage::sparse::state_type state) const;


};
}  // namespace storage
}  // namespace storm

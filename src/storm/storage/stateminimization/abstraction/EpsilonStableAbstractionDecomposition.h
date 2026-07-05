#pragma once

#include <deque>
#include <fstream>
#include <random>
#include <utility>

#include "storm/adapters/IntervalForward.h"
#include "storm/exceptions/InvalidOptionException.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/settings/modules/AbstractionSettings.h"
#include "storm/storage/Distribution.h"
#include "storm/storage/stateminimization/BaseDecomposition.h"
#include "storm/utility/interval.h"

namespace storm {
namespace storage {
namespace abstraction {

template<typename ModelType>
class EpsilonStableAbstractionDecomposition : public BaseDecomposition<ModelType> {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef typename ModelType::RewardModelType RewardModelType;

    struct EpsilonStableAbstractionOptions : BaseDecomposition<ModelType>::BaseOptions {
       public:
        EpsilonStableAbstractionOptions() = default;

        EpsilonStableAbstractionOptions(double epsilon);

        /*!
         * Creates an object representing the options necessary to obtain the quotient while still preserving
         * the given formula.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formula The formula that is to be preserved.
         */
        EpsilonStableAbstractionOptions(ModelType const& model, storm::logic::Formula const& formula, double epsilon);

        /*!
         * Creates an object representing the options necessary to obtain the smallest quotient while still
         * preserving the given formulas.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formulas The formulas that need to be preserved.
         */
        EpsilonStableAbstractionOptions(ModelType const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas, double epsilon);

        double getEpsilon() const {
            return this->epsilon;
        }

        void setEpsilon(double epsilon) {
            this->epsilon = epsilon;
        }

       private:
        /// Represents the epsilon allowed for the abstraction.
        double epsilon;
    };

    struct StateGroup {
       public:
        StateGroup() = default;

        StateGroup(std::vector<storm::storage::sparse::state_type> groupStates,
                   std::span<storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const> groupDistributions)
            : states(std::move(groupStates)), choiceDistributions(groupDistributions.begin(), groupDistributions.end()) {}

        std::vector<storm::storage::sparse::state_type> const& getStates() const {
            return states;
        }

        void addState(storm::storage::sparse::state_type state) {
            states.push_back(state);
        }

        std::vector<storage::Distribution<ValueType, storm::storage::sparse::state_type>> const& getDistributions() const {
            return choiceDistributions;
        }

        storage::Distribution<ValueType, storm::storage::sparse::state_type> const& getDistributionAtOffset(std::uint_fast64_t choiceOffset) const {
            return choiceDistributions[choiceOffset];
        }

        void setDistributionAtOffset(std::uint_fast64_t choiceOffset,
                                     storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> distribution) {
            choiceDistributions[choiceOffset] = std::move(distribution);
        }

        std::size_t getNumberOfChoices() {
            return choiceDistributions.size();
        }

       private:
        /// States of the group.
        std::vector<storm::storage::sparse::state_type> states;

        /// One enhanced group distribution per enabled choice/action position.
        /// For deterministic models, this vector has size 1.
        std::vector<storage::Distribution<ValueType, storm::storage::sparse::state_type>> choiceDistributions;
    };

    /*!
     * Decomposes the given model into equivalence classes satisfying the epsilon-stability.
     *
     * @param model The model to decompose.
     * @param options The options to use during for the decomposition.
     */
    EpsilonStableAbstractionDecomposition(ModelType const& model, EpsilonStableAbstractionOptions const& options);

   protected:
    void computeInitialPartition() override;

    void refinePartition() override;

    bool shouldBuildQuotient() const override;

    void buildQuotientFromPartition() override;

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

   private:
    void splitInitialPartitionBasedOnLabelsAndRewards();

    void splitInitialPartitionBasedOnRewards();

    void splitInitialPartitionBasedOnRewards(std::vector<ValueType> const& rewardVector);

    void splitInitialPartitionBasedOnActionRewards(std::vector<std::set<ValueType>> const& rewardVector);

    void splitInitialPartitionBasedOnActionSets();

    void initializeChoiceDistributions();

    void recomputeChoiceDistribution(uint_fast64_t choice);

    void recomputeChoiceDistributionsOfState(storm::storage::sparse::state_type state);

    std::span<storm::storage::Distribution<ValueType, storm::storage::sparse::state_type> const> getChoiceDistributionsOfState(
        storm::storage::sparse::state_type state) const;

    std::pair<std::uint_fast64_t, std::uint_fast64_t> getChoiceRangeOfState(storm::storage::sparse::state_type state) const;

    ValueType computeFeasibleIntervalBasedOnAggregatedIntervals(ValueType intervalToSplitter, ValueType intervalToOtherBlocks) const;

    /*!
     * Performs the epsilon-stable abstraction on the model and thereby merges suitable states which satisfy
     * the epsilon-stability criterion. If required, the quotient model is built and may be retrieved using
     * getQuotient().
     * @param epsilon
     */
    void performEpsilonStableAbstractionRefinement(double epsilon);

    void refineBlockBasedOnEpsilonSignature(std::span<uint64_t const> subBlock, std::deque<typename stateminimization::Partition::Block>& blocksQueue,
                                            stateminimization::Partition::BlockSet& enqueuedBlocks, double epsilon);

    EpsilonStableAbstractionOptions options;

    /// One distribution per state for IMCs/IDTMCs, one distribution for every choice/action per state for IMDPs.
    std::vector<storm::storage::Distribution<ValueType, storm::storage::sparse::state_type>> originalChoiceDistributions;
};

}  // namespace abstraction
}  // namespace storage
}  // namespace storm

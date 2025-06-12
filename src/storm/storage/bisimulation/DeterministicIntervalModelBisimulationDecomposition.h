#ifndef STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H
#define STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H

#include <storm/storage/bisimulation/Partition.h>
#include "storm/storage/bisimulation/BisimulationDecomposition.h"
#include "storm/models/sparse/DeterministicModel.h"
#include <storm/storage/geometry/Polytope.h>

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

   private:
    void postProcessInitialPartition();

    std::shared_ptr<storm::storage::geometry::Polytope<storm::RationalNumber>> create2DPolytope(storm::RationalNumber c1LowerBound, storm::RationalNumber c1UpperBound,
                     storm::RationalNumber c2LowerBound, storm::RationalNumber c2UpperBound);
};

}

}



#endif  // STORM_DETERMINISTICINTERVALMODELBISIMULATIONDECOMPOSITION_H

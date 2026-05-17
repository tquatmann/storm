#pragma once

#include <deque>
#include <random>

#include "storm/exceptions/AbortException.h"
#include "storm/exceptions/InvalidOptionException.h"
#include "storm/models/sparse/Model.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/BisimulationSettings.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/solver/OptimizationDirection.h"
#include "storm/storage/Decomposition.h"
#include "storm/storage/StateBlock.h"
#include "storm/storage/sparse/StateType.h"
#include "storm/storage/stateminimization/BaseDecomposition.h"
#include "storm/storage/stateminimization/bisimulation/BisimulationType.h"
#include "storm/storage/stateminimization/bisimulation/Signature.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/constants.h"

namespace storm {
namespace logic {
class Formula;
}

namespace storage {

// TODO: Do we even need this class?
// TODO: As there is not much implementation still left here, we could just let the implementations of the bisimulation directly implement the BaseDecomposition
// TODO: class.
/*!
 * This class is the superclass of all decompositions of a sparse model into its bisimulation quotient.
 */
template<typename ModelType>
class BisimulationDecomposition : public BaseDecomposition<ModelType> {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef typename ModelType::RewardModelType RewardModelType;

    // A class that offers the possibility to customize the bisimulation.
    struct BisimulationOptions : public BaseDecomposition<ModelType>::BaseOptions {
        // Creates an object representing the default values for all options.
        BisimulationOptions();

        /*!
         * Creates an object representing the options necessary to obtain the quotient while still preserving
         * the given formula.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formula The formula that is to be preserved.
         */
        BisimulationOptions(ModelType const& model, storm::logic::Formula const& formula);

        /*!
         * Creates an object representing the options necessary to obtain the smallest quotient while still
         * preserving the given formulas.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formulas The formulas that need to be preserved.
         */
        BisimulationOptions(ModelType const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);

        /**
         * Sets the bisimulation type. If the bisimulation type is set to weak,
         * we also change the bounded flag (as bounded properties are not preserved under
         * weak bisimulation).
         */
        void setType(BisimulationType t) {
            if (t == BisimulationType::Weak) {
                STORM_LOG_WARN_COND(!this->isBounded(), "Weak bisimulation does not preserve bounded properties.");
                STORM_LOG_WARN_COND(!this->isDiscounted(), "Weak bisimulation does not preserve discounted properties.");
                this->setBounded(false);
                this->setDiscounted(false);
            }
            type = t;
        }

        BisimulationType getType() const {
            return this->type;
        }

       private:
        /// A flag that indicates whether a strong or a weak bisimulation is to be computed.
        BisimulationType type;
    };

    /*!
     * Decomposes the given model into equivalence classes of a bisimulation.
     *
     * @param model The model to decompose.
     * @param options The options to use during for the decomposition.
     */
    BisimulationDecomposition(ModelType const& model, BisimulationOptions const& options);

    virtual ~BisimulationDecomposition() = default;

   protected:
    /*!
     * Decomposes the given model into equivalence classes of a bisimulation.
     *
     * @param model The model to decompose.
     * @param backwardTransitions The backward transitions of the model.
     * @param options The options to use during for the decomposition.
     */
    BisimulationDecomposition(ModelType const& model, storm::storage::SparseMatrix<ValueType> const& backwardTransitions, BisimulationOptions const& options);

    void computeInitialPartition() override;

    void refinePartition() override;

    bool shouldBuildQuotient() const override;

    /*!
     * Performs the partition refinement on the model and thereby computes the equivalence classes under strong
     * bisimulation equivalence. If required, the quotient model is built and may be retrieved using
     * getQuotient().
     */
    void performPartitionRefinement();

    /*!
     * Performs the signature refinement on the model and thereby computes the equivalence classes under strong
     * bisimulation equivalence. If required, the quotient model is built and may be retrieved using
     * getQuotient().
     */
    void performSignatureRefinement();

    /*!
     * Computes the signature of the given state with respect to the given partition.
     * @param state input state whose signature shall be computed.
     * @param currentPartition current currentPartition to compute the signature.
     * @return hash value of the state's signature.
     */
    storm::storage::bisimulation::Signature<typename ModelType::ValueType> computeStateSignature(
        storm::storage::sparse::state_type state, storm::storage::stateminimization::Partition const& currentPartition) const;

    /*!
     * Computes a hash value based of the signature of the given state with respect to the given partition.
     * @param state input state whose signature shall be computed.
     * @return hash value of the state's signature.
     */
    std::size_t computeStateSignatureHash(storm::storage::sparse::state_type state) const;

    /*!
     * Refines the partition by considering the given splitter. All blocks that become potential splitters
     * because of this refinement, are marked as splitters and inserted into the splitter vector.
     *
     * @param splitterBlock The splitter to use.
     */
    virtual void refinePartitionBasedOnSplitter(std::span<uint64_t const> splitterBlock,
                                                std::deque<typename stateminimization::Partition::Block>& splitterQueue,
                                                stateminimization::Partition::BlockSet& enqueuedSplitterBlocks) = 0;

    /*!
     * Initializes the initial partition based on all respected labels.
     */
    void initializeLabelBasedPartition();

    /*!
     * Creates the measure-driven initial partition for reaching psi states from phi states.
     */
    void initializeMeasureDrivenPartition();

    /*!
     * A function that can initialize auxiliary data structures. It is called after initializing the initial partition.
     */
    virtual void initialize();

    void splitInitialPartitionBasedOnActionSets();

    /*!
     * Computes the set of states with probability 0/1 for satisfying phi until psi. This is used for the measure
     * driven initial partition.
     *
     * @return The states with probability 0 and 1.
     */
    virtual std::pair<storm::storage::BitVector, storm::storage::BitVector> getStatesWithProbability01() = 0;

    /*!
     * Splits the initial partition based on the (unique) reward model of the current model.
     */
    virtual void splitInitialPartitionBasedOnRewards();

    /*!
     * Splits the initial partition based on the given reward vector.
     */
    virtual void splitInitialPartitionBasedOnRewards(std::vector<ValueType> const& rewardVector);

    /*!
     * Splits the initial partition based on the given vector of action rewards.
     */
    virtual void splitInitialPartitionBasedOnActionRewards(std::vector<std::set<ValueType>> const& rewardVector);

    /*!
     * Constructs the blocks of the decomposition object based on the current partition.
     */
    void extractDecompositionBlocks();

    BisimulationOptions options;
};

inline BisimulationType resolveBisimulationTypeChoice(BisimulationTypeChoice c) {
    switch (c) {
        case BisimulationTypeChoice::Strong:
            return BisimulationType::Strong;
        case BisimulationTypeChoice::Weak:
            return BisimulationType::Weak;
        case BisimulationTypeChoice::FromSettings:
            if (storm::settings::getModule<storm::settings::modules::BisimulationSettings>().isWeakBisimulationSet()) {
                return BisimulationType::Weak;
            } else {
                return BisimulationType::Strong;
            }
    }
    return BisimulationType::Strong;
}

}  // namespace storage
}  // namespace storm

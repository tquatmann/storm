#pragma once

#include "storm/adapters/IntervalForward.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/logic/FormulaInformation.h"
#include "storm/logic/Formulas.h"
#include "storm/logic/FragmentSpecification.h"
#include "storm/modelchecker/propositional/SparsePropositionalModelChecker.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/Model.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/CoreSettings.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/stateminimization/Partition.h"
#include "storm/utility/constants.h"

namespace storm::storage {

template<typename ModelType>
class BaseDecomposition {
   public:
    typedef typename ModelType::ValueType ValueType;
    typedef typename ModelType::RewardModelType RewardModelType;

    BaseDecomposition(ModelType const& model, storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::IntervalBaseType<ValueType> tolerance);

    // A class that offers the possibility to customize the decomposition.
    struct BaseOptions {
       public:
        // Creates an object representing the default values for all options.
        BaseOptions();

        /*!
         * Creates an object representing the options necessary to obtain the quotient while still preserving
         * the given formula.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formula The formula that is to be preserved.
         */
        BaseOptions(ModelType const& model, storm::logic::Formula const& formula);

        /*!
         * Creates an object representing the options necessary to obtain the smallest quotient while still
         * preserving the given formulas.
         *
         * @param model The model for which the quotient model shall be computed. This needs to be given in order to
         * derive a suitable initial partition.
         * @param formulas The formulas that need to be preserved.
         */
        BaseOptions(ModelType const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas);

        void preserveFormula(storm::logic::Formula const& formula);

        void preserveSingleFormula(ModelType const& model, storm::logic::Formula const& formula);

        void checkAndSetMeasureDrivenInitialPartition(ModelType const& model, storm::logic::Formula const& formula);

        void addToRespectedAtomicPropositions(std::vector<std::shared_ptr<storm::logic::AtomicExpressionFormula const>> const& expressions,
                                              std::vector<std::shared_ptr<storm::logic::AtomicLabelFormula const>> const& labels);

        bool getKeepRewards() const {
            return this->keepRewards;
        }

        void setKeepRewards(bool keepRewards) {
            this->keepRewards = keepRewards;
        }

        bool isOptimizationDirectionSet() const {
            return static_cast<bool>(optimalityType);
        }

        OptimizationDirection getOptimizationDirection() const {
            STORM_LOG_ASSERT(optimalityType, "Optimality type not set.");
            return optimalityType.value();
        }

        storm::IntervalBaseType<ValueType> getTolerance() const {
            return storm::NumberTraits<storm::IntervalBaseType<ValueType>>::IsExact
                       ? storm::utility::zero<storm::IntervalBaseType<ValueType>>()
                       : storm::utility::convertNumber<storm::IntervalBaseType<ValueType>>(
                             storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision());
        }

        bool isBounded() const {
            return bounded;
        }

        bool isDiscounted() const {
            return discounted;
        }

        void setBounded(bool value) {
            bounded = value;
        }

        void setDiscounted(bool value) {
            discounted = value;
        }

        /// A flag that indicates whether a measure driven initial partition is to be used. If this flag is set
        /// to true, the two optional pairs phiStatesAndLabel and psiStatesAndLabel must be set. Then, the
        /// measure driven initial partition wrt. to the states phi and psi is taken.
        bool measureDrivenInitialPartition;
        std::optional<storm::storage::BitVector> phiStates;
        std::optional<storm::storage::BitVector> psiStates;

        /// An optional set of strings that indicate which of the atomic propositions of the model are to be
        /// respected and which may be ignored. If not given, all atomic propositions of the model are respected.
        std::optional<std::set<std::string>> respectedAtomicPropositions;

        /// A flag that governs whether the quotient model is actually built or only the decomposition is computed.
        bool buildQuotient;

       private:
        std::optional<OptimizationDirection> optimalityType;

        /// A flag that indicates whether the state-rewards of the model are to be respected (and should
        /// be kept in the quotient model, if one is built).
        bool keepRewards;

        /// A flag that indicates whether step-bounded properties are to be preserved. This may only be set to true
        /// when computing strong bisimulation equivalence.
        bool bounded;

        /// A flag that indicates whether discounted properties are to be preserved. This may only be set to true
        /// when computing strong bisimulation equivalence.
        bool discounted;
    };

    /*!
     * Computes the decomposition of the model into equivalence classes. If requested, a quotient model is built.
     */
    void computeDecomposition();

    /*!
     * Retrieves the quotient of the model under the computed decomposition.
     *
     * @return The quotient model.
     */
    std::shared_ptr<ModelType> getQuotient() const;

   protected:
    using Clock = std::chrono::high_resolution_clock;
    using Duration = Clock::duration;

    virtual void computeInitialPartition() = 0;

    virtual void refinePartition() = 0;

    virtual bool shouldBuildQuotient() const = 0;

    virtual void buildQuotientFromPartition() = 0;

    template<typename Callback>
    typename Clock::duration measure(Callback&& callback) {
        auto start = Clock::now();
        std::forward<Callback>(callback)();
        return Clock::now() - start;
    }

    void printStatistics(typename Clock::duration initialPartitionTime, typename Clock::duration refinementTime, typename Clock::duration quotientBuildTime,
                         typename Clock::duration totalTime) const;

    /// The model to decompose.
    ModelType const& model;

    /// The backward transitions of the model.
    storm::storage::SparseMatrix<ValueType> backwardTransitions;

    /// The current partition.
    storm::storage::stateminimization::Partition partition;

    /// The quotient, if it was build. Otherwise, a null pointer.
    std::shared_ptr<ModelType> quotient;

    /// Map of representative states of absorbing blocks. A single entry represents: <first state of block, representative state>
    std::map<uint64_t, uint64_t> absorbingBlocks;

    /// A comparator used for comparing the distances of constants.
    storm::utility::ConstantsComparator<storm::IntervalBaseType<ValueType>> comparator;
};

}  // namespace storm::storage

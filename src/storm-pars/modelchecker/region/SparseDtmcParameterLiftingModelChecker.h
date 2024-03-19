#pragma once

#include <memory>
#include <optional>
#include <vector>

#include "storm-pars/modelchecker/instantiation/SparseDtmcInstantiationModelChecker.h"
#include "storm-pars/modelchecker/region/SparseParameterLiftingModelChecker.h"
#include "storm-pars/transformer/ParameterLifter.h"
#include "storm/logic/FragmentSpecification.h"
#include "storm/solver/MinMaxLinearEquationSolver.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/Scheduler.h"

namespace storm {
namespace modelchecker {

template<typename ParametricType, typename ConstantType>
class OrderBasedMonotonicityBackend;

template<typename SparseModelType, typename ConstantType>
class SparseDtmcParameterLiftingModelChecker : public SparseParameterLiftingModelChecker<SparseModelType, ConstantType> {
   public:
    using ParametricType = typename SparseModelType::ValueType;
    using CoefficientType = typename RegionModelChecker<ParametricType>::CoefficientType;
    using VariableType = typename RegionModelChecker<ParametricType>::VariableType;
    using Valuation = typename RegionModelChecker<ParametricType>::Valuation;

    SparseDtmcParameterLiftingModelChecker();
    SparseDtmcParameterLiftingModelChecker(std::unique_ptr<storm::solver::MinMaxLinearEquationSolverFactory<ConstantType>>&& solverFactory);
    virtual ~SparseDtmcParameterLiftingModelChecker() = default;

    virtual bool canHandle(std::shared_ptr<storm::models::ModelBase> parametricModel,
                           CheckTask<storm::logic::Formula, ParametricType> const& checkTask) const override;

    virtual void specify(Environment const& env, std::shared_ptr<storm::models::ModelBase> parametricModel,
                         CheckTask<storm::logic::Formula, ParametricType> const& checkTask,
                         std::optional<RegionSplitEstimateKind> generateRegionSplitEstimates = std::nullopt,
                         std::shared_ptr<MonotonicityBackend<ParametricType>> monotonicityBackend = {}, bool allowModelSimplifications = true) override;

    std::optional<storm::storage::Scheduler<ConstantType>> getCurrentMinScheduler();
    std::optional<storm::storage::Scheduler<ConstantType>> getCurrentMaxScheduler();

    virtual bool isRegionSplitEstimateKindSupported(RegionSplitEstimateKind kind,
                                                    CheckTask<storm::logic::Formula, ParametricType> const& checkTask) const override;
    virtual RegionSplitEstimateKind getDefaultRegionSplitEstimateKind(CheckTask<storm::logic::Formula, ParametricType> const& checkTask) const override;
    virtual std::vector<CoefficientType> obtainRegionSplitEstimates(std::set<VariableType> const& relevantParameters) const override;

    virtual bool isMonotonicitySupported(MonotonicityBackend<ParametricType> const& backend,
                                         CheckTask<storm::logic::Formula, ParametricType> const& checkTask) const override;

   protected:
    virtual void specifyBoundedUntilFormula(const CheckTask<storm::logic::BoundedUntilFormula, ConstantType>& checkTask) override;
    virtual void specifyUntilFormula(Environment const& env, CheckTask<storm::logic::UntilFormula, ConstantType> const& checkTask) override;
    virtual void specifyReachabilityRewardFormula(Environment const& env, CheckTask<storm::logic::EventuallyFormula, ConstantType> const& checkTask) override;
    virtual void specifyCumulativeRewardFormula(const CheckTask<storm::logic::CumulativeRewardFormula, ConstantType>& checkTask) override;

    virtual storm::modelchecker::SparseInstantiationModelChecker<SparseModelType, ConstantType>& getInstantiationChecker() override;
    virtual storm::modelchecker::SparseInstantiationModelChecker<SparseModelType, ConstantType>& getInstantiationCheckerSAT() override;
    virtual storm::modelchecker::SparseInstantiationModelChecker<SparseModelType, ConstantType>& getInstantiationCheckerVIO() override;

    virtual std::unique_ptr<CheckResult> computeQuantitativeValues(Environment const& env, AnnotatedRegion<ParametricType>& region,
                                                                   storm::solver::OptimizationDirection const& dirForParameters) override;

    void computeStateValueDeltaRegionSplitEstimates(std::vector<ConstantType> const& quantitativeResult, std::vector<uint64_t> const& schedulerChoices,
                                                    storm::storage::ParameterRegion<ParametricType> const& region,
                                                    storm::solver::OptimizationDirection const& dirForParameters);

    virtual void reset() override;

   private:
    void computeSchedulerDeltaSplitEstimates(std::vector<ConstantType> const& quantitativeResult, std::vector<uint64_t> const& schedulerChoices,
                                             storm::storage::ParameterRegion<ParametricType> const& region,
                                             storm::solver::OptimizationDirection const& dirForParameters);

    bool isOrderBasedMonotonicityBackend() const;
    OrderBasedMonotonicityBackend<ParametricType, ConstantType>& getOrderBasedMonotonicityBackend();

    bool isValueDeltaRegionSplitEstimates() const;

    storm::storage::BitVector maybeStates;
    std::vector<ConstantType> resultsForNonMaybeStates;
    std::optional<uint64_t> stepBound;

    std::unique_ptr<storm::modelchecker::SparseDtmcInstantiationModelChecker<SparseModelType, ConstantType>> instantiationChecker;
    std::unique_ptr<storm::modelchecker::SparseDtmcInstantiationModelChecker<SparseModelType, ConstantType>> instantiationCheckerSAT;
    std::unique_ptr<storm::modelchecker::SparseDtmcInstantiationModelChecker<SparseModelType, ConstantType>> instantiationCheckerVIO;

    std::unique_ptr<storm::transformer::ParameterLifter<ParametricType, ConstantType>> parameterLifter;
    std::unique_ptr<storm::solver::MinMaxLinearEquationSolverFactory<ConstantType>> solverFactory;
    bool solvingRequiresUpperRewardBounds;

    // Results from the most recent solver call.
    std::optional<std::vector<uint64_t>> minSchedChoices, maxSchedChoices;
    std::vector<ConstantType> x;
    std::optional<ConstantType> lowerResultBound, upperResultBound;

    std::map<VariableType, ConstantType> cachedRegionSplitEstimates;
};
}  // namespace modelchecker
}  // namespace storm

#pragma once

#include <memory>

#include "storm/modelchecker/parametric/RegionCheckResult.h"
#include "storm/modelchecker/parametric/SparseInstantiationModelChecker.h"
#include "storm/modelchecker/parametric/SparseParameterLiftingModelChecker.h"
#include "storm/modelchecker/CheckTask.h"
#include "storm/storage/ParameterRegion.h"
#include "storm/utility/Stopwatch.h"

namespace storm {
    namespace modelchecker{
        namespace parametric{
            
            struct ParameterLiftingSettings {
                ParameterLiftingSettings();
                
                bool applyExactValidation;
            };

            template<typename SparseModelType, typename ConstantType>
            class ParameterLifting {

            public:
                
                typedef typename storm::storage::ParameterRegion<typename SparseModelType::ValueType>::CoefficientType CoefficientType;

                ParameterLifting(SparseModelType const& parametricModel);

                ParameterLiftingSettings const& getSettings() const;
                void setSettings(ParameterLiftingSettings const& newSettings);
                
                void specifyFormula(CheckTask<storm::logic::Formula, typename SparseModelType::ValueType> const& checkTask);
                
                /*!
                 * Analyzes the given region by means of parameter lifting.
                 * We first check whether there is one point in the region for which the property is satisfied/violated.
                 * If the given initialResults already indicates that there is such a point, this step is skipped.
                 * Then, we check whether ALL points in the region violate/satisfy the property
                 *
                 */
                RegionCheckResult analyzeRegion(storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, RegionCheckResult const& initialResult = RegionCheckResult::Unknown, bool sampleVerticesOfRegion = false);
                /*!
                 * Similar to analyze region but additionaly invokes exact parameter lifting to validate results AllSat or AllViolated
                 */
                RegionCheckResult analyzeRegionExactValidation(storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, RegionCheckResult const& initialResult = RegionCheckResult::Unknown);
                
                /*!
                 * Iteratively refines the region until parameter lifting yields a conclusive result (AllSat or AllViolated).
                 * The refinement stops as soon as the fraction of the area of the subregions with inconclusive result is less then the given threshold
                 */
                std::vector<std::pair<storm::storage::ParameterRegion<typename SparseModelType::ValueType>, RegionCheckResult>> performRegionRefinement(storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& region, CoefficientType const& threshold);
                
                static std::string visualizeResult(std::vector<std::pair<storm::storage::ParameterRegion<typename SparseModelType::ValueType>, RegionCheckResult>> const& result, storm::storage::ParameterRegion<typename SparseModelType::ValueType> const& parameterSpace, typename storm::storage::ParameterRegion<typename SparseModelType::ValueType>::VariableType const& x, typename storm::storage::ParameterRegion<typename SparseModelType::ValueType>::VariableType const& y);
                
            protected:
                SparseModelType const& getConsideredParametricModel() const;
                
                virtual void initializeUnderlyingCheckers() = 0;
                virtual void simplifyParametricModel(CheckTask<storm::logic::Formula, typename SparseModelType::ValueType> const& checkTask) = 0;
                virtual void applyHintsToExactChecker() = 0;
                
                SparseModelType const& parametricModel;
                ParameterLiftingSettings settings;
                std::unique_ptr<CheckTask<storm::logic::Formula, typename SparseModelType::ValueType>> currentCheckTask;
                std::shared_ptr<storm::logic::Formula const> currentFormula;
                std::shared_ptr<SparseModelType> simplifiedModel;
                

                std::unique_ptr<SparseParameterLiftingModelChecker<SparseModelType, ConstantType>> parameterLiftingChecker;
                std::unique_ptr<SparseParameterLiftingModelChecker<SparseModelType, storm::RationalNumber>> exactParameterLiftingChecker; // todo: use template argument instead of rational number
                //std::unique_ptr<SparseParameterLiftingModelChecker<SparseModelType, CoefficientType>> exactParameterLiftingChecker;
                std::unique_ptr<SparseInstantiationModelChecker<SparseModelType, ConstantType>> instantiationChecker;
                
                mutable storm::utility::Stopwatch initializationStopwatch, instantiationCheckerStopwatch, parameterLiftingCheckerStopwatch;
            };
    
        } //namespace parametric
    } //namespace modelchecker
} //namespace storm

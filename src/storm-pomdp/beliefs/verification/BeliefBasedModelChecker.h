#pragma once

#include "storm-pomdp/beliefs/verification/PropertyInformation.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelCheckerOptions.h"
#include "storm-pomdp/storage/BeliefExplorationBounds.h"

namespace storm {
class Environment;

namespace pomdp::beliefs {

template<typename PomdpModelType, typename BeliefValueType = typename PomdpModelType::ValueType,
         typename BeliefMdpValueType = typename PomdpModelType::ValueType>
class BeliefBasedModelChecker {
   public:
    explicit BeliefBasedModelChecker(PomdpModelType const& pomdp);

    BeliefMdpValueType checkUnfold(storm::Environment const& env, PropertyInformation const& propertyInformation,
                                   storm::pomdp::storage::PreprocessingPomdpValueBounds<BeliefMdpValueType> const& valueBounds = {});

    BeliefMdpValueType checkDiscretize(storm::Environment const& env, PropertyInformation const& propertyInformation, uint64_t resolution, bool useDynamic,
                                       storm::pomdp::storage::PreprocessingPomdpValueBounds<BeliefMdpValueType> const& valueBounds = {});

   private:
    PomdpModelType const& inputPomdp;
};
}  // namespace pomdp::beliefs
}  // namespace storm
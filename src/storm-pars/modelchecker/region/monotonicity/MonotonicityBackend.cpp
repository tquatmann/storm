#include "storm-pars/modelchecker/region/monotonicity/MonotonicityBackend.h"

#include "storm/utility/macros.h"

namespace storm::modelchecker {

template<typename ParametricType>
void MonotonicityBackend<ParametricType>::setMonotoneParameter(VariableType const& parameter, MonotonicityKind const& kind) {
    STORM_LOG_ASSERT(storm::analysis::isMonotone(kind), "Monotonicity kind must be either increasing, decreasing or constant.");
    globallyKnownMonotonicityInformation[parameter] = kind;
}

template<typename ParametricType>
void MonotonicityBackend<ParametricType>::initializeMonotonicity(AnnotatedRegion<ParametricType>& region) {
    typename AnnotatedRegion<ParametricType>::DefaultMonotonicityAnnotation annotation;
    annotation.globalMonotonicity = std::make_shared<storm::analysis::MonotonicityResult<VariableType>>();
    bool allMonotone = true;
    for (auto const& parameter : region.region.getVariables()) {
        if (auto findRes = globallyKnownMonotonicityInformation.find(parameter); findRes != globallyKnownMonotonicityInformation.end()) {
            annotation.globalMonotonicity->addMonotonicityResult(parameter, findRes->second);
            annotation.globalMonotonicity->setDoneForVar(parameter);
        } else {
            allMonotone = false;
        }
    }
    annotation.globalMonotonicity->setAllMonotonicity(allMonotone);
    annotation.globalMonotonicity->setDone();
    region.monotonicityAnnotation = annotation;
}

template<typename ParametricType>
void MonotonicityBackend<ParametricType>::updateMonotonicity(AnnotatedRegion<ParametricType>&) {
    // Nothing to do here
}

template<typename ParametricType>
bool MonotonicityBackend<ParametricType>::requiresInteractionWithRegionModelChecker() const {
    return false;
}

template<typename ParametricType>
bool MonotonicityBackend<ParametricType>::recommendModelSimplifications() const {
    return true;
}

template<typename ParametricType>
std::map<typename MonotonicityBackend<ParametricType>::VariableType, typename MonotonicityBackend<ParametricType>::MonotonicityKind>
MonotonicityBackend<ParametricType>::getOptimisticMonotonicityApproximation(AnnotatedRegion<ParametricType> const& region) {
    auto annotatedResult = region.getDefaultMonotonicityAnnotation();
    STORM_LOG_ASSERT(annotatedResult.has_value() && annotatedResult->globalMonotonicity, "Default monotonicity annotation with a result must be present.");
    std::map<VariableType, MonotonicityKind> result;
    for (auto const& parameter : region.region.getVariables()) {
        if (auto monRes = annotatedResult->globalMonotonicity->getMonotonicity(parameter); storm::analysis::isMonotone(monRes)) {
            result.emplace(parameter, monRes);
        }
    }
    return result;
}

template class MonotonicityBackend<storm::RationalFunction>;

}  // namespace storm::modelchecker
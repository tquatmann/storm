#include "storm-pomdp/beliefs/verification/BeliefBasedModelChecker.h"

#include <memory>

#include "storm-pomdp/beliefs/abstraction/FreudenthalTriangulationBeliefAbstraction.h"
#include "storm-pomdp/beliefs/exploration/BeliefExploration.h"
#include "storm-pomdp/beliefs/exploration/BeliefMdpBuilder.h"
#include "storm-pomdp/beliefs/storage/Belief.h"
#include "storm-pomdp/beliefs/utility/types.h"

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/api/verification.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/models/sparse/Pomdp.h"
#include "storm/utility/OptionalRef.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

template<typename PomdpModelType, typename BeliefValueType, typename BeliefMdpValueType>

BeliefBasedModelChecker<PomdpModelType, BeliefValueType, BeliefMdpValueType>::BeliefBasedModelChecker(PomdpModelType const& pomdp) : inputPomdp(pomdp) {
    STORM_LOG_ERROR_COND(inputPomdp.isCanonic(), "Input Pomdp is not known to be canonic. This might lead to unexpected verification results.");
}

template<typename PomdpModelType, typename BeliefType, typename BeliefMdpValueType>
BeliefMdpValueType checkUnfoldOrDiscretize(storm::Environment const& env, PomdpModelType const& pomdp, PropertyInformation const& propertyInformation,
                                           storm::pomdp::storage::PreprocessingPomdpValueBounds<BeliefMdpValueType> const& valueBounds,
                                           storm::OptionalRef<FreudenthalTriangulationBeliefAbstraction<BeliefType>> abstraction = {}) {
    STORM_LOG_ASSERT(propertyInformation.kind == PropertyInformation::Kind::ReachabilityProbability ||
                         propertyInformation.kind == PropertyInformation::Kind::ExpectedTotalReachabilityReward,
                     "Unexpected kind of property.");

    STORM_PRINT_AND_LOG("Constructing the belief MDP...\n");

    // First, explore the beliefs and its successors
    using BeliefExplorationType = BeliefExploration<BeliefMdpValueType, PomdpModelType, BeliefType>;
    storm::utility::Stopwatch swExplore(true);
    BeliefExplorationType exploration(pomdp);
    auto info = exploration.initializeStandardExploration();
    // TODO: implement a TerminationCallback (e.g. based on the number of explored states) that prevents the exploration from running indefinitely in case of
    // infinite belief MDPs
    // typename BeliefExplorationType::TerminationCallback terminationCallback = [&info]() {};
    if (propertyInformation.kind == PropertyInformation::Kind::ExpectedTotalReachabilityReward) {
        typename BeliefExplorationType::TerminalBeliefCallback terminalBeliefCallback =
            [&propertyInformation](BeliefType const& belief) -> std::optional<BeliefMdpValueType> {
            if (propertyInformation.targetObservations.count(belief.observation()) > 0) {
                return storm::utility::zero<BeliefMdpValueType>();
            } else {
                return std::nullopt;
            };
        };
        exploration.resumeExploration(info, terminalBeliefCallback, {}, propertyInformation.rewardModelName.value(), abstraction);
    } else {
        typename BeliefExplorationType::TerminalBeliefCallback terminalBeliefCallback =
            [&propertyInformation](BeliefType const& belief) -> std::optional<BeliefMdpValueType> {
            if (propertyInformation.targetObservations.count(belief.observation()) > 0) {
                return storm::utility::one<BeliefMdpValueType>();
            } else {
                return std::nullopt;
            };
        };
        exploration.resumeExploration(info, terminalBeliefCallback, {}, storm::NullRef, abstraction);
    }
    swExplore.stop();
    STORM_LOG_WARN_COND(!info.queue.hasNext(), "Exploration stopped before all states were explored.");

    // Second, build the Belief MDP from the exploration information
    storm::utility::Stopwatch swBuild(true);
    std::function<BeliefMdpValueType(BeliefType const&)> computeCutOffValue = [&valueBounds](BeliefType const& belief) {
        // TODO: use value bounds to compute the cut-off value
        assert(false);
        return storm::utility::zero<BeliefMdpValueType>();
    };
    auto beliefMdp = buildBeliefMdp(info, propertyInformation, computeCutOffValue);
    swBuild.stop();
    beliefMdp->printModelInformationToStream(std::cout);

    // Finally, perform model checking on the belief MDP.
    storm::utility::Stopwatch swCheck(true);
    auto formula = createFormulaForBeliefMdp(propertyInformation);
    storm::modelchecker::CheckTask<storm::logic::Formula, BeliefMdpValueType> task(*formula, true);
    std::unique_ptr<storm::modelchecker::CheckResult> res(storm::api::verifyWithSparseEngine<BeliefMdpValueType>(env, beliefMdp, task));
    swCheck.stop();
    STORM_PRINT_AND_LOG("Time for exploring beliefs: " << swExplore << ".\n");
    STORM_PRINT_AND_LOG("Time for building the belief MDP: " << swBuild << ".\n");
    STORM_PRINT_AND_LOG("Time for analyzing the belief MDP: " << swCheck << ".\n");
    STORM_LOG_ASSERT(res, "Model checking of belief MDP did not return any result.");
    STORM_LOG_ASSERT(res->isExplicitQuantitativeCheckResult(), "Model checking of belief MDP did not return result of expectedd type.");
    STORM_LOG_ASSERT(beliefMdp->getInitialStates().getNumberOfSetBits() == 1, "Unexpected number of initial states for belief Mdp.");
    auto const initState = beliefMdp->getInitialStates().getNextSetIndex(0);
    return res->asExplicitQuantitativeCheckResult<BeliefMdpValueType>()[initState];
}

template<typename PomdpModelType, typename BeliefValueType, typename BeliefMdpValueType>
BeliefMdpValueType BeliefBasedModelChecker<PomdpModelType, BeliefValueType, BeliefMdpValueType>::checkUnfold(
    storm::Environment const& env, PropertyInformation const& propertyInformation,
    storm::pomdp::storage::PreprocessingPomdpValueBounds<BeliefMdpValueType> const& valueBounds) {
    return checkUnfoldOrDiscretize<PomdpModelType, Belief<BeliefValueType>, BeliefMdpValueType>(env, inputPomdp, propertyInformation, valueBounds);
}

template<typename PomdpModelType, typename BeliefValueType, typename BeliefMdpValueType>
BeliefMdpValueType BeliefBasedModelChecker<PomdpModelType, BeliefValueType, BeliefMdpValueType>::checkDiscretize(
    storm::Environment const& env, PropertyInformation const& propertyInformation, uint64_t resolution, bool useDynamic,
    storm::pomdp::storage::PreprocessingPomdpValueBounds<BeliefMdpValueType> const& valueBounds) {
    std::vector<BeliefValueType> observationResolutionVector(inputPomdp.getNrObservations(), storm::utility::convertNumber<BeliefValueType>(resolution));
    auto mode = useDynamic ? FreudenthalTriangulationMode::Dynamic : FreudenthalTriangulationMode::Static;
    FreudenthalTriangulationBeliefAbstraction<Belief<BeliefValueType>> abstraction(observationResolutionVector, mode);
    return checkUnfoldOrDiscretize<PomdpModelType, Belief<BeliefValueType>, BeliefMdpValueType>(env, inputPomdp, propertyInformation, valueBounds, abstraction);
}

// TODO: Check which instantiations are actually necessary / reasonable.
template class BeliefBasedModelChecker<storm::models::sparse::Pomdp<double>, double, double>;
template class BeliefBasedModelChecker<storm::models::sparse::Pomdp<double>, storm::RationalNumber, double>;
template class BeliefBasedModelChecker<storm::models::sparse::Pomdp<storm::RationalNumber>, double, storm::RationalNumber>;
template class BeliefBasedModelChecker<storm::models::sparse::Pomdp<storm::RationalNumber>, storm::RationalNumber, storm::RationalNumber>;
}  // namespace storm::pomdp::beliefs

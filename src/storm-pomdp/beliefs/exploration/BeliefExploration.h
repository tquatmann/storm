#pragma once

#include <functional>
#include <optional>
#include <set>

#include "storm-pomdp/beliefs/exploration/ExplorationInformation.h"

#include "storm-pomdp/beliefs/exploration/FirstStateNextStateGenerator.h"
#include "storm-pomdp/beliefs/utility/types.h"

#include "storm/utility/OptionalRef.h"
#include "storm/utility/vector.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
class FreudenthalTriangulationBeliefAbstraction;

template<typename BeliefMdpValueType, typename PomdpType, typename BeliefType>
class BeliefExploration {
   public:
    using TerminationCallback = std::function<bool()>;
    using TerminalBeliefCallback = std::function<std::optional<BeliefMdpValueType>(BeliefType const&)>;

    BeliefExploration(PomdpType const& pomdp);

    using StandardExplorationInformation = ExplorationInformation<BeliefMdpValueType, BeliefType>;
    StandardExplorationInformation initializeStandardExploration(ExplorationQueueOrder const explorationQueueOrder = ExplorationQueueOrder::Unordered);
    void resumeExploration(StandardExplorationInformation& info, TerminalBeliefCallback const& terminalBeliefCallback = {},
                           TerminationCallback const& terminationCallback = {}, storm::OptionalRef<std::string const> rewardModelName = {},
                           storm::OptionalRef<FreudenthalTriangulationBeliefAbstraction<BeliefType>> abstraction = {});

   private:
    template<typename InfoType, typename NextStateHandleType>
    bool performExploration(InfoType& info, NextStateHandleType&& exploreNextStates, TerminalBeliefCallback const& terminalBeliefCallback,
                            TerminationCallback const& terminationCallback);

    storm::pomdp::beliefs::FirstStateNextStateGenerator<PomdpType, BeliefType> firstStateNextStateGenerator;
};
}  // namespace storm::pomdp::beliefs
#pragma once

#include <set>

#include "storm-pomdp/beliefs/exploration/BeliefExplorationMode.h"

#include "storm-pomdp/beliefs/abstraction/RewardBoundedBeliefSplitter.h"
#include "storm-pomdp/beliefs/exploration/FirstStateNextStateGenerator.h"
#include "storm-pomdp/beliefs/storage/BeliefCollector.h"
#include "storm-pomdp/beliefs/utility/types.h"
#include "storm/utility/OptionalRef.h"
#include "storm/utility/vector.h"

namespace storm::pomdp::builder {

namespace detail {

template<typename ValueType>
struct StandardExplorationTransition {
    ValueType probability;
    storm::pomdp::beliefs::BeliefId targetBelief;
};

template<typename ValueType>
struct MultiRewardAugmentedExplorationTransition {
    ValueType probability;
    storm::pomdp::beliefs::BeliefId targetBelief;
    std::vector<ValueType> reward;
};
}  // namespace detail

template<typename ValueType, BeliefExplorationMode Mode>
using BeliefExplorationTransition = std::conditional_t<Mode == BeliefExplorationMode::Standard, detail::StandardExplorationTransition<ValueType>,
                                                       detail::MultiRewardAugmentedExplorationTransition<ValueType>>;

template<typename ValueType, BeliefExplorationMode Mode>
class BeliefExplorationMatrix {
   public:
    BeliefExplorationMatrix() {
        rowIndications.push_back(0u);
        rowGroupIndices.push_back(0u);
    }
    std::size_t rows() const {
        return rowIndications.size() - 1;
    }
    std::size_t groups() const {
        return rowGroupIndices.size() - 1;
    }

    void endCurrentRow() {
        rowIndications.push_back(transitions.size());
    };

    void endCurrentRowGroup() {
        rowGroupIndices.push_back(rowIndications.size() - 1);
    };

    std::vector<BeliefExplorationTransition<ValueType, Mode>> transitions;
    std::vector<uint64_t> rowIndications;
    std::vector<uint64_t> rowGroupIndices;
};

}  // namespace storm::pomdp::builder
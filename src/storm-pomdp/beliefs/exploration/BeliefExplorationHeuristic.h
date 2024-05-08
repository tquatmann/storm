#pragma once

#include <set>

#include "storm/storage/SparseMatrix.h"

#include "storm-pomdp/beliefs/storage/BeliefCollector.h"
#include "storm-pomdp/beliefs/utility/types.h"

namespace storm::pomdp::builder {

// TODO: rename to  backend
template<typename BeliefType>
class BeliefExplorationHeuristic {
   public:
    using BeliefId = storm::pomdp::beliefs::BeliefId;
    bool hasNext() const {
        return !queue.empty();
    }

    bool discover(BeliefId const& id, BeliefType const& belief) {
        if (isTerminalBelief(belief)) {
            return false;
        }
        return queue.insert(id).second;
    }

    BeliefId popNext() {
        assert(hasNext());
        BeliefId id = *queue.begin();
        queue.erase(queue.begin());
        return id;
    }

    void setTerminalObservations(std::set<storm::pomdp::beliefs::BeliefObservationType> const& observations) {
        terminalObservations = observations;
    }

   private:
    bool isTerminalBelief(BeliefType const& belief) {
        return terminalObservations.count(belief.observation()) != 0;
    }

    std::set<storm::pomdp::beliefs::BeliefId> queue;
    std::set<storm::pomdp::beliefs::BeliefObservationType> terminalObservations;
};

}  // namespace storm::pomdp::builder
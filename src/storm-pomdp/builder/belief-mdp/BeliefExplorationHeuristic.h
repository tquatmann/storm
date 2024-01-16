#pragma once

#include <set>

#include "storm/storage/SparseMatrix.h"

#include "storm-pomdp/storage/beliefs/BeliefCollector.h"
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"

namespace storm::pomdp::builder {

class BeliefExplorationHeuristic {
   public:
    using BeliefId = storm::pomdp::beliefs::BeliefId;
    bool hasNext() const {
        return !queue.empty();
    }

    bool discover(BeliefId const& id) {
        return queue.insert(id).second;
    }

    BeliefId popNext() {
        assert(hasNext());
        BeliefId id = *queue.begin();
        queue.erase(queue.begin());
        return id;
    }

   private:
    std::set<storm::pomdp::beliefs::BeliefId> queue;
};

}  // namespace storm::pomdp::builder
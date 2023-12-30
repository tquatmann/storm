#pragma once
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"

#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
class BeliefCollector {
   public:
    /*!
     * Indicates that this collector gives IDs in a consecutive way, ensuring that all Ids are less than `getNumberOfBeliefIds()`
     */
    bool HasConsecutiveIds = true;

    bool isEqual(BeliefId const& first, BeliefId const& second) const {
        return first == second || getBelief(first) == getBelief(second);
    }

    BeliefId getNumberOfBeliefIds() const {
        return gatheredBeliefs.size();
    }

    BeliefType const& getBelief(BeliefId const& id) const {
        STORM_LOG_ASSERT(id < gatheredBeliefs.size(), "Unexpected belief id " << id << ". Ids are in [0," << getNumberOfBeliefIds() << ").");
        return gatheredBeliefs[id];
    }

    /*!
     * @return the ID of the given belief. Assumes that the belief is present in the collector.
     */
    BeliefId getId(BeliefType const& belief) const {
        STORM_LOG_ASSERT(belief.observation() < beliefToIdMap.size(),
                         "Unknown belief observation " << belief.observation() << ". Obervations are in [0," << beliefToIdMap.size());
        STORM_LOG_ASSERT(hasBelief(belief), "Belief " << belief.toString() << " is not present in this collector.");
        return beliefToIdMap[belief.observation()].at(belief);
    }

    bool hasBelief(BeliefType const& belief) const {
        return belief.observation() < beliefToIdMap.size() && beliefToIdMap[belief.observation()].count(belief) > 0;
    }

    /*!
     * @return The Id of the given belief if it is present in the collector. `InvalidBeliefId` otherwise.
     */
    BeliefId getIdOptional(BeliefType const& belief) const {
        if (auto const& obs = belief.observation(); obs < beliefToIdMap.size()) {
            if (auto findRes = beliefToIdMap[obs].find(belief); findRes != beliefToIdMap[obs].end()) {
                return findRes->second;
            }
        }
        return InvalidBeliefId;
    }

    /*!
     * Checks if the given belief has already been collected before.
     *  - If yes, the Id is returned.
     *  - If not, the belief is collected and a fresh Id is allocated.
     */
    BeliefId getOrAddBeliefId(BeliefType&& belief) {
        if (auto id = getIdOptional(belief); id != InvalidBeliefId) {
            return id;
        }
        return addBeliefId(std::move(belief));
    }

    /*!
     * Adds a fresh belief and returns its ID. Assumes that this belief is not already present.
     * @return
     */
    BeliefId addBeliefId(BeliefType&& inputBelief) {
        auto const id = gatheredBeliefs.size();
        gatheredBeliefs.push_back(std::move(inputBelief));
        auto const& belief = gatheredBeliefs.back();
        auto const& obs = belief.observation();
        if (obs >= beliefToIdMap.size()) {
            beliefToIdMap.resize(obs + 1);
        }
        beliefToIdMap[obs].emplace(belief, id);
        return id;
    }

   private:
    std::vector<BeliefType> gatheredBeliefs;
    std::vector<std::unordered_map<BeliefType, BeliefId, typename BeliefType::BeliefHash>> beliefToIdMap;
};
}  // namespace storm::pomdp::beliefs
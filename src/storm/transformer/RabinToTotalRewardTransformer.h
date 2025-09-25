#pragma once

#include <memory>
#include <string>
#include <vector>

#include "storm/transformer/LTLToRabinObjectiveTransformer.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotImplementedException.h"

namespace storm {

namespace storage {
class BitVector;
}

namespace transformer {

/*!
 * Demerges MECs
 */
class RabinToTotalRewardTransformer {
   public:
    using RabinCondition = LTLToRabinObjectiveTransformer::RabinCondition;
    using RabinObjective = LTLToRabinObjectiveTransformer::RabinObjective;

    template<typename SparseModelType>
    struct ReturnType {
        std::shared_ptr<SparseModelType> model;
    };

    /*!
     * Produces a model with total reward objectives where every MEC of the input model is replaced by a set of new MECs and a proxy state that decides which
     * new MEC to enter (or whether to move on to another component). Each such new MEC is actually a sub-EC of the original MEC.
     * Entering a new MEC gives a reward of 1 for each quantitativeRabinObjective that is satisfied when each state is visited infinitely often.
     * The i'th total reward objective corresponds to the i'th quantitativeRabinObjective.
     * All new MECs satisfy the almostSureRabinObjectives (assuming that each state in that MEC is visited infinitely often).
     * States that can not satisfy all allmostSureRabinObjectives and choices leading to such states are removed.
     * If that includes all initial states, a nullptr model is returned.
     *
     * @return
     */
    template<typename SparseModelType>
    static ReturnType<SparseModelType> transform(SparseModelType const& model, std::vector<RabinObjective> const& quantitativeRabinObjectives,
                                                 std::vector<std::string> const& quantitativeTotalRewardModelNames,
                                                 std::vector<RabinObjective> const& almostSureRabinObjectives) {
        STORM_LOG_ASSERT(quantitativeRabinObjectives.size() == quantitativeTotalRewardModelNames.size(),
                         "Inconsistent number of quantitative rabin objectives.");

        STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Rabin to total reward transformation is not implemented for sparse models yet.");
    }
};
}  // namespace transformer
}  // namespace storm
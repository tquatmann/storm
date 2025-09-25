#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "storm/logic/FormulasForwardDeclarations.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotImplementedException.h"
namespace storm {

namespace storage {
class BitVector;
}

namespace transformer {

/*!
 * Incorporates Memory into the state space of the given model, that is
 * the resulting model is the crossproduct of of the given model plus
 * some type of memory structure
 */
class LTLToRabinObjectiveTransformer {
   public:
    struct RabinCondition {
        std::vector<std::string> fin, inf;  // state labels that should be visited finitely / infinitely often (read as conjunction)
    };
    using RabinObjective = std::vector<RabinCondition>;  // disjunction of Rabin conditions

    template<typename SparseModelType>
    struct ReturnType {
        std::shared_ptr<SparseModelType> model;
        std::vector<RabinObjective> rabinObjectives;
    };

    /*!
     * Takes as input a model and a set of LTL formulas and produces a model with equivalent Rabin acceptance conditions for each formula.
     * @param env
     * @param model
     * @param ltlFormulas
     * @return
     */
    template<typename SparseModelType>
    static ReturnType<SparseModelType> transform(SparseModelType const& model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
                                                 std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& subformulaCallback) {
        STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Rabin to total reward transformation is not implemented for sparse models yet.");

        /*
         template<class SparseModelType>
std::shared_ptr<SparseModelType> incorporateLTLMemory(Environment const& env, SparseModelType const& model,
                                             std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) {
if constexpr (std::is_same<SparseModelType, models::sparse::Mdp<typename SparseModelType::ValueType>>()) {
std::vector<std::shared_ptr<storm::logic::Formula const>> ltlFormulas;
std::vector<std::shared_ptr<storm::logic::PathFormula const>> test;

for (auto const& subFormula : formulas) {
   STORM_LOG_THROW(subFormula->isOperatorFormula(), storm::exceptions::NotSupportedException,
                   "The given Formula " << *subFormula << " is not supported.");

   if (subFormula->isProbabilityOperatorFormula()) {
       ltlFormulas.push_back(subFormula);
   }
}

} else {
throw storm::exceptions::NotImplementedException() << "The LTL memory incorporation is only implemented for MDPs.";
}
}

        */
    }
};
}  // namespace transformer
}  // namespace storm
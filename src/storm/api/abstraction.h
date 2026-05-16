#pragma once

#include "storm/logic/Formula.h"
#include "storm/models/sparse/Ctmc.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/storage/stateminimization/abstraction/AbstractionType.h"
#include "storm/storage/stateminimization/abstraction/EpsilonStableAbstractionDecomposition.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/NotSupportedException.h"

namespace storm {
namespace api {

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> performAbstractionMinimization(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> model, std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
    storm::storage::stateminimization::abstraction::AbstractionType abstractionType, bool graphPreserving = true, double epsilon = 0.0) {
    STORM_LOG_THROW((model->isOfType(storm::models::ModelType::Dtmc) || model->isOfType(storm::models::ModelType::Mdp)) && storm::IsIntervalType<ValueType>,
                    storm::exceptions::NotSupportedException, "Epsilon-stable abstraction minimization is currently only available for IDTMCs and IMDPs.");

    if (model->isOfType(storm::models::ModelType::Dtmc)) {
        auto options = typename storm::storage::abstraction::EpsilonStableAbstractionDecomposition<
            storm::models::sparse::Dtmc<ValueType>>::EpsilonStableAbstractionOptions(*model->template as<storm::models::sparse::Dtmc<ValueType>>(), formulas,
                                                                                     epsilon);

        storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Dtmc<ValueType>> abstractionDecomposition(
            *model->template as<storm::models::sparse::Dtmc<ValueType>>(), options);
        abstractionDecomposition.computeDecomposition();
        return abstractionDecomposition.getQuotient();
    } else {
        auto options =
            typename storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<ValueType>>::EpsilonStableAbstractionOptions(
                *model->template as<storm::models::sparse::Mdp<ValueType>>(), formulas, epsilon);

        storm::storage::abstraction::EpsilonStableAbstractionDecomposition<storm::models::sparse::Mdp<ValueType>> abstractionDecomposition(
            *model->template as<storm::models::sparse::Mdp<ValueType>>(), options);
        abstractionDecomposition.computeDecomposition();
        return abstractionDecomposition.getQuotient();
    }
}

}  // namespace api
}  // namespace storm
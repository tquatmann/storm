#pragma once

#include <numeric>
#include "storm/automata/AcceptanceCondition.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

namespace storm::automata {
using BooleanHoaExpression = cpphoafparser::BooleanExpression<cpphoafparser::AtomAcceptance>;

namespace detail {

/*!
 * Adds the given offset to all acceptance set indices in the given expression.
 */
inline BooleanHoaExpression::ptr addOffsetAcceptanceSets(BooleanHoaExpression::ptr const expr, uint64_t offset) {
    if (expr->isTRUE() || expr->isFALSE()) {
        return expr;
    }

    if (expr->isAtom()) {
        auto const& atom = expr->getAtom();
        auto const offsetAtom = std::make_shared<cpphoafparser::AtomAcceptance>(atom.getType(), atom.getAcceptanceSet() + offset, atom.isNegated());
        return BooleanHoaExpression::Atom(offsetAtom);
    }

    if (expr->isAND() || expr->isOR()) {
        auto left = addOffsetAcceptanceSets(expr->getLeft(), offset);
        auto right = addOffsetAcceptanceSets(expr->getRight(), offset);
        return expr->isAND() ? left & right : left | right;
    }

    if (expr->isNOT()) {
        auto left = addOffsetAcceptanceSets(expr->getLeft(), offset);
        return !left;
    }

    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unhandled expression when adding offset to acceptance sets");
}
}  // namespace detail

/*!
 * Calls the given callback for each combination of the acceptance conditions provided in the constructor
 * (i.e., for each subset of the acceptance conditions).
 * The callback is called with two parameters:
 * 1. The combined acceptance condition (with acceptance sets merged and the acceptance expression being the conjunction of the enabled conditions)
 * 2. A BitVector indicating which of the original acceptance conditions are enabled (true) and which are disabled (false)
 *
 * @param acceptanceConditions
 * @param callback
 */
inline void forEachAcceptanceCombination(std::vector<AcceptanceCondition::ptr> const& acceptanceConditions,
                                         std::function<void(AcceptanceCondition::ptr, storm::storage::BitVector const&)> const& callback) {
    STORM_LOG_ASSERT(acceptanceConditions.size() > 0, "No acceptance conditions provided");
    std::vector<BooleanHoaExpression::ptr> acceptanceExpressions;
    std::vector<storm::storage::BitVector> acceptanceSets;
    for (const auto& acceptanceCondition : acceptanceConditions) {
        acceptanceExpressions.push_back(detail::addOffsetAcceptanceSets(acceptanceCondition->getAcceptanceExpression(), acceptanceSets.size()));
        for (uint64_t i = 0; i < acceptanceCondition->getNumberOfAcceptanceSets(); ++i) {
            acceptanceSets.push_back(acceptanceCondition->getAcceptanceSet(i));
        }
    }

    auto accCond = std::make_shared<AcceptanceCondition>(std::move(acceptanceSets), BooleanHoaExpression::True());
    storm::storage::BitVector enabledConditions(acceptanceConditions.size(), false);
    STORM_LOG_INFO("Calling callback with acceptance expression: " << accCond->getAcceptanceExpression() << "\n\tand enabled conditions " << enabledConditions);
    callback(accCond, enabledConditions);  // first call with all conditions disabled
    do {
        enabledConditions.increment();
        bool first = true;
        BooleanHoaExpression::ptr accExprPtr;
        for (auto condIndex : enabledConditions) {
            if (first) {
                accExprPtr = acceptanceExpressions[condIndex];
                first = false;
            } else {
                accExprPtr = accExprPtr & acceptanceExpressions[condIndex];
            }
        }
        accCond->setAcceptanceExpression(accExprPtr);
        accCond->convertToDNF();
        STORM_LOG_INFO("Calling callback with acceptance expression: " << accCond->getAcceptanceExpression() << "\n\tand enabled conditions "
                                                                       << enabledConditions);
        callback(accCond, enabledConditions);
    } while (!enabledConditions.full());
}

}  // namespace storm::automata

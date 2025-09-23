#pragma once

#include <numeric>
#include "storm/automata/AcceptanceCondition.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

namespace storm::automata {

namespace detail {

using HoaBoolExpr = cpphoafparser::BooleanExpression<cpphoafparser::AtomAcceptance>;

HoaBoolExpr::ptr toNegationNormalForm(HoaBoolExpr::ptr expr) {
    if (expr->isTRUE() || expr->isFALSE() || expr->isAtom()) {
        return expr;
    }
    if (expr->isAND()) {
        return toNegationNormalForm(expr->getLeft()) & toNegationNormalForm(expr->getRight());
    }
    if (expr->isOR()) {
        return toNegationNormalForm(expr->getLeft()) | toNegationNormalForm(expr->getRight());
    }
    if (expr->isNOT()) {
        auto subExpr = expr->getLeft();
        if (subExpr->isTRUE())
            return HoaBoolExpr::False();
        if (subExpr->isFALSE())
            return HoaBoolExpr::True();
        if (subExpr->isAtom()) {
            using AtomAcceptance = cpphoafparser::AtomAcceptance;
            auto const& atom = subExpr->getAtom();
            auto negatedAtomType = atom.getType() == AtomAcceptance::TEMPORAL_FIN ? AtomAcceptance::TEMPORAL_INF : AtomAcceptance::TEMPORAL_FIN;
            auto const negatedAtom = std::make_shared<AtomAcceptance>(negatedAtomType, atom.getAcceptanceSet(), atom.isNegated());
            return HoaBoolExpr::Atom(negatedAtom);
        }
        if (subExpr->isNOT()) {
            return toNegationNormalForm(subExpr->getLeft());
        }
        if (subExpr->isAND()) {
            // De Morgan: !(A & B) -> !A | !B
            return toNegationNormalForm(!subExpr->getLeft() | !subExpr->getRight());
        }
        if (subExpr->isOR()) {
            // De Morgan: !(A | B) -> !A & !B
            return toNegationNormalForm(!subExpr->getLeft() & !subExpr->getRight());
        }
    }
    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unhandled expression in negation normal form conversion");
}

HoaBoolExpr::ptr fromNNFToDNF(HoaBoolExpr::ptr expr, bool topLevel = true) {
    if (expr->isAND()) {
        auto left = fromNNFToDNF(expr->getLeft());
        auto right = fromNNFToDNF(expr->getRight());

        // TODO: this can be optimized, e.g., check for clauses that are always false (x & !x) or that subsume each other (x & y & z | x & y)
        if (left->isOR()) {
            // Distributive law: (A | B) & C -> (A & C) | (B & C)
            return fromNNFToDNF((left->getLeft() & right) | (left->getRight() & right));
        }
        if (right->isOR()) {
            // Distributive law: A & (B | C) -> (A & B) | (A & C)
            return fromNNFToDNF((right->getLeft() & left) | (right->getRight() & left));
        }
        // Nothing to distribute, just combine
        return left & right;
    }

    if (expr->isOR()) {
        return fromNNFToDNF(expr->getLeft()) | fromNNFToDNF(expr->getRight());
    }

    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unhandled expression in disjunctive normal form conversion");
}

HoaBoolExpr::ptr toDisjunctiveNormalForm(HoaBoolExpr::ptr expr) {
    return fromNNFToDNF(toNegationNormalForm(expr));
}

/*!
 * Adds the given offset to all acceptance set indices in the given expression.
 */
HoaBoolExpr::ptr addOffsetAcceptanceSets(HoaBoolExpr::ptr const expr, uint64_t offset) {
    if (expr->isTRUE() || expr->isFALSE()) {
        return expr;
    }

    if (expr->isAtom()) {
        auto const& atom = expr->getAtom();
        auto const offsetAtom = std::make_shared<cpphoafparser::AtomAcceptance>(atom.getType(), atom.getAcceptanceSet() + offset, atom.isNegated());
        return HoaBoolExpr::Atom(offsetAtom);
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
void forEachAcceptanceCombination(std::vector<AcceptanceCondition::ptr> const& acceptanceConditions,
                                  std::function<void(AcceptanceCondition::ptr, storm::storage::BitVector const&)> const& callback) {
    STORM_LOG_ASSERT(acceptanceConditions.size() > 0, "No acceptance conditions provided");
    std::vector<detail::HoaBoolExpr::ptr> acceptanceExpressions;
    std::vector<storm::storage::BitVector> acceptanceSets;
    for (const auto& acceptanceCondition : acceptanceConditions) {
        acceptanceExpressions.push_back(detail::addOffsetAcceptanceSets(acceptanceCondition->getAcceptanceExpression(), acceptanceSets.size()));
        for (uint64_t i = 0; i < acceptanceCondition->getNumberOfAcceptanceSets(); ++i) {
            acceptanceSets.push_back(acceptanceCondition->getAcceptanceSet(i));
        }
    }

    detail::HoaBoolExpr::ptr accExprPtr = detail::HoaBoolExpr::True();
    auto accCond = std::make_shared<AcceptanceCondition>(std::move(acceptanceSets), accExprPtr);
    storm::storage::BitVector enabledConditions(acceptanceConditions.size(), false);
    STORM_LOG_INFO("Calling callback with acceptance expression: " << accExprPtr->toString() << "\n\tand enabled conditions " << enabledConditions);
    callback(accCond, enabledConditions);  // first call with all conditions disabled
    do {
        enabledConditions.increment();
        bool first = true;
        for (auto condIndex : enabledConditions) {
            if (first) {
                accExprPtr = acceptanceExpressions[condIndex];
                first = false;
            } else {
                accExprPtr = accExprPtr & acceptanceExpressions[condIndex];
            }
        }
        accExprPtr = detail::toDisjunctiveNormalForm(accExprPtr);
        STORM_LOG_INFO("Calling callback with acceptance expression: " << accExprPtr->toString() << "\n\tand enabled conditions " << enabledConditions);
        callback(accCond, enabledConditions);
    } while (!enabledConditions.full());
}

}  // namespace storm::automata

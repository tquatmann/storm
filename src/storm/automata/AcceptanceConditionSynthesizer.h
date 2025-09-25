#pragma once

#include <numeric>
#include "storm/automata/AcceptanceCondition.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

namespace storm::automata {
using BooleanHoaExpression = cpphoafparser::BooleanExpression<cpphoafparser::AtomAcceptance>;

namespace detail {

inline BooleanHoaExpression::ptr toNegationNormalForm(BooleanHoaExpression::ptr expr) {
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
            return BooleanHoaExpression::False();
        if (subExpr->isFALSE())
            return BooleanHoaExpression::True();
        if (subExpr->isAtom()) {
            using AtomAcceptance = cpphoafparser::AtomAcceptance;
            auto const& atom = subExpr->getAtom();
            auto negatedAtomType = atom.getType() == AtomAcceptance::TEMPORAL_FIN ? AtomAcceptance::TEMPORAL_INF : AtomAcceptance::TEMPORAL_FIN;
            auto const negatedAtom = std::make_shared<AtomAcceptance>(negatedAtomType, atom.getAcceptanceSet(), atom.isNegated());
            return BooleanHoaExpression::Atom(negatedAtom);
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

inline BooleanHoaExpression::ptr fromNNFToDNF(BooleanHoaExpression::ptr expr, bool topLevel = true) {
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

inline bool isConjunctionOfLiterals(BooleanHoaExpression::ptr expr) {
    // Catch (negated) literals
    auto e = expr->isNOT() ? expr->getLeft() : expr;
    if (e->isTRUE() || e->isFALSE() || e->isAtom()) {
        return true;
    }

    if (expr->isAND()) {
        return isConjunctionOfLiterals(expr->getLeft()) && isConjunctionOfLiterals(expr->getRight());
    }
    return false;
}

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

inline bool isDisjunctiveNormalForm(BooleanHoaExpression::ptr expr) {
    if (expr->isOR()) {
        return isDisjunctiveNormalForm(expr->getLeft()) && isDisjunctiveNormalForm(expr->getRight());
    } else {
        return detail::isConjunctionOfLiterals(expr);
    }
}

inline BooleanHoaExpression::ptr toDisjunctiveNormalForm(BooleanHoaExpression::ptr expr) {
    return detail::fromNNFToDNF(detail::toNegationNormalForm(expr));
}

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
        accCond->setAcceptanceExpression(toDisjunctiveNormalForm(accExprPtr));
        STORM_LOG_INFO("Calling callback with acceptance expression: " << accCond->getAcceptanceExpression() << "\n\tand enabled conditions "
                                                                       << enabledConditions);
        callback(accCond, enabledConditions);
    } while (!enabledConditions.full());
}

}  // namespace storm::automata

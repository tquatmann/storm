#include "AcceptanceCondition.h"

#include "storm/exceptions/InvalidOperationException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

namespace storm {
namespace automata {

AcceptanceCondition::AcceptanceCondition(std::size_t numberOfStates, unsigned int numberOfAcceptanceSets, acceptance_expr::ptr acceptance)
    : acceptanceSets(numberOfStates, storm::storage::BitVector(numberOfStates)), acceptance(acceptance) {
    // intentionally left empty
}

AcceptanceCondition::AcceptanceCondition(std::vector<storm::storage::BitVector>&& acceptanceSets, acceptance_expr::ptr acceptance)
    : acceptanceSets(std::move(acceptanceSets)), acceptance(acceptance) {
    // intentionally left empty
}

unsigned int AcceptanceCondition::getNumberOfAcceptanceSets() const {
    return acceptanceSets.size();
}

storm::storage::BitVector& AcceptanceCondition::getAcceptanceSet(unsigned int index) {
    return acceptanceSets.at(index);
}

const storm::storage::BitVector& AcceptanceCondition::getAcceptanceSet(unsigned int index) const {
    return acceptanceSets.at(index);
}

AcceptanceCondition::acceptance_expr::ptr AcceptanceCondition::getAcceptanceExpression() const {
    return acceptance;
}

void AcceptanceCondition::setAcceptanceExpression(acceptance_expr::ptr expr) {
    acceptance = expr;
}

bool AcceptanceCondition::isAccepting(const storm::storage::StateBlock& scc) const {
    return isAccepting(scc, acceptance);
}

bool AcceptanceCondition::isAccepting(const storm::storage::StateBlock& scc, acceptance_expr::ptr expr) const {
    switch (expr->getType()) {
        case acceptance_expr::EXP_AND:
            return isAccepting(scc, expr->getLeft()) && isAccepting(scc, expr->getRight());
        case acceptance_expr::EXP_OR:
            return isAccepting(scc, expr->getLeft()) || isAccepting(scc, expr->getRight());
        case acceptance_expr::EXP_NOT:
            return !isAccepting(scc, expr->getLeft());
        case acceptance_expr::EXP_TRUE:
            return true;
        case acceptance_expr::EXP_FALSE:
            return false;
        case acceptance_expr::EXP_ATOM: {
            const cpphoafparser::AtomAcceptance& atom = expr->getAtom();
            const storm::storage::BitVector& acceptanceSet = acceptanceSets.at(atom.getAcceptanceSet());
            bool negated = atom.isNegated();
            bool rv;
            switch (atom.getType()) {
                case cpphoafparser::AtomAcceptance::TEMPORAL_INF:
                    rv = false;
                    for (auto& state : scc) {
                        if (acceptanceSet.get(state)) {
                            rv = true;
                            break;
                        }
                    }
                    break;
                case cpphoafparser::AtomAcceptance::TEMPORAL_FIN:
                    rv = true;
                    for (auto& state : scc) {
                        if (acceptanceSet.get(state)) {
                            rv = false;
                            break;
                        }
                    }
                    break;
            }

            return (negated ? !rv : rv);
        }
    }

    throw std::runtime_error("Missing case statement");
}

namespace detail {
using BooleanHoaExpression = cpphoafparser::BooleanExpression<cpphoafparser::AtomAcceptance>;

void extractFromDNFRecursion(AcceptanceCondition::acceptance_expr::ptr e, std::vector<std::vector<BooleanHoaExpression::ptr>>& dnf, bool topLevel) {
    if (topLevel) {
        if (e->isOR()) {
            if (e->getLeft()->isOR()) {
                extractFromDNFRecursion(e->getLeft(), dnf, true);
            } else {
                dnf.emplace_back();
                extractFromDNFRecursion(e->getLeft(), dnf, false);
            }

            if (e->getRight()->isOR()) {
                extractFromDNFRecursion(e->getRight(), dnf, true);
            } else {
                dnf.emplace_back();
                extractFromDNFRecursion(e->getRight(), dnf, false);
            }
        } else {
            dnf.emplace_back();
            extractFromDNFRecursion(e, dnf, false);
        }
    } else {
        if (e->isOR() || e->isNOT()) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Acceptance condition is not in DNF");
        } else if (e->isAND()) {
            extractFromDNFRecursion(e->getLeft(), dnf, false);
            extractFromDNFRecursion(e->getRight(), dnf, false);
        } else {
            dnf.back().push_back(e);
        }
    }
}

BooleanHoaExpression::ptr toNegationNormalForm(BooleanHoaExpression::ptr expr) {
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

BooleanHoaExpression::ptr fromNNFToDNF(BooleanHoaExpression::ptr expr, bool topLevel = true) {
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

bool isDNF(BooleanHoaExpression::ptr expr, bool topLevel = true) {
    if (topLevel && expr->isOR()) {  // OR only allowed at top level
        return isDNF(expr->getLeft()) && isDNF(expr->getRight());
    }
    // Reaching this point means we are no longer at the top level
    if (expr->isAND()) {
        return isDNF(expr->getLeft(), false) && isDNF(expr->getRight(), false);
    }

    // Catch (negated) literals
    auto e = expr->isNOT() ? expr->getLeft() : expr;
    return e->isAtom() || e->isTRUE() || e->isFALSE();
}

}  // namespace detail

std::vector<std::vector<AcceptanceCondition::acceptance_expr::ptr>> AcceptanceCondition::extractFromDNF() const {
    std::vector<std::vector<AcceptanceCondition::acceptance_expr::ptr>> dnf;

    detail::extractFromDNFRecursion(getAcceptanceExpression(), dnf, true);

    return dnf;
}

bool AcceptanceCondition::isInDNF() const {
    return detail::isDNF(getAcceptanceExpression());
}

void AcceptanceCondition::convertToDNF() {
    setAcceptanceExpression(detail::fromNNFToDNF(detail::toNegationNormalForm(getAcceptanceExpression())));
}

AcceptanceCondition::ptr AcceptanceCondition::lift(std::size_t productNumberOfStates, std::function<std::size_t(std::size_t)> mapping) const {
    AcceptanceCondition::ptr lifted(new AcceptanceCondition(productNumberOfStates, getNumberOfAcceptanceSets(), acceptance));
    for (uint64_t i = 0; i < getNumberOfAcceptanceSets(); ++i) {
        const storm::storage::BitVector& set = getAcceptanceSet(i);
        storm::storage::BitVector& liftedSet = lifted->getAcceptanceSet(i);

        for (std::size_t prodState = 0; prodState < productNumberOfStates; prodState++) {
            if (set.get(mapping(prodState))) {
                liftedSet.set(prodState);
            }
        }
    }

    return lifted;
}

}  // namespace automata
}  // namespace storm

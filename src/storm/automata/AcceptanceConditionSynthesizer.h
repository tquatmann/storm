#pragma once

#include <numeric>

#include "../utility/logging.h"
#include "AcceptanceCondition.h"

namespace cpphoafparser {

BooleanExpression<AtomAcceptance>::ptr toDNF(
    BooleanExpression<AtomAcceptance>::ptr expr) {
    using ptr = BooleanExpression<AtomAcceptance>::ptr;

    if (expr->isTRUE() || expr->isFALSE() || expr->isAtom()) {
        return expr;
    }

    if (expr->isNOT()) {
        // Negationen behandeln und durch De Morgan'sche Regeln vereinfachen
        ptr subExpr = expr->getLeft();
        if (subExpr->isTRUE()) return BooleanExpression<AtomAcceptance>::False();
        if (subExpr->isFALSE()) return BooleanExpression<AtomAcceptance>::True();
        if (subExpr->isAtom()) {
            auto const& atom = subExpr->getAtom();
            auto negatedAtomType = atom.getType() == AtomAcceptance::TEMPORAL_FIN ?
                AtomAcceptance::TEMPORAL_INF :
                AtomAcceptance::TEMPORAL_FIN;
            auto const negatedAtom = std::make_shared<AtomAcceptance>(negatedAtomType, atom.getAcceptanceSet(), atom.isNegated());
            return BooleanExpression<AtomAcceptance>::Atom(negatedAtom);
        }
        if (subExpr->isAND()) {
            // De Morgan: !(A & B) -> !A | !B
            return toDNF(!subExpr->getLeft() | !subExpr->getRight());
        }
        if (subExpr->isOR()) {
            // De Morgan: !(A | B) -> !A & !B
            return toDNF(!subExpr->getLeft() & !subExpr->getRight());
        }
    }

    if (expr->isAND()) {
        ptr left = toDNF(expr->getLeft());
        ptr right = toDNF(expr->getRight());

        if (left->isOR()) {
            // Distributivgesetz: (A | B) & C -> (A & C) | (B & C)
            return toDNF((left->getLeft() & right) | (left->getRight() & right));
        }
        if (right->isOR()) {
            // Distributivgesetz: A & (B | C) -> (A & B) | (A & C)
            return toDNF((right->getLeft() & left) | (right->getRight() & left));
        }
        // Keine Änderung erforderlich
        return left & right;
    }

    if (expr->isOR()) {
        // Rekursive Anwendung auf Disjunktionen
        return toDNF(expr->getLeft()) | toDNF(expr->getRight());
    }

    throw std::logic_error("Unbekannter Operator im Ausdruck");
}

BooleanExpression<AtomAcceptance>::ptr addOffsetAcceptanceSets(BooleanExpression<AtomAcceptance>::ptr const expr, uint64_t offset) {
    if (expr->isTRUE() || expr->isFALSE()) {
        return expr;
    }

    if (expr->isAtom()) {
        auto const& atom = expr->getAtom();
        auto const offsetAtom = std::make_shared<AtomAcceptance>(atom.getType(), atom.getAcceptanceSet() + offset, atom.isNegated());
        return BooleanExpression<AtomAcceptance>::Atom(offsetAtom);
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

    throw std::logic_error("Unbekannter Operator im Ausdruck");
}

}  // namespace cpphoafparser


namespace storm {
namespace automata {
class AcceptanceConditionSynthesizer {
public:
    AcceptanceConditionSynthesizer(std::vector<AcceptanceCondition::ptr> acceptanceConditions): acceptanceConditions(acceptanceConditions) {}

 /**
  * Makes the numbering of acceptance sets unique, e.g., [Inf(0), Inf(0)] -> [Inf(0), Inf(1)]
  *
  * @param acceptanceConditions The vector of acceptance conditions for the LTL objectives
  */
 static std::vector<AcceptanceCondition::ptr> liftAcceptanceSets(std::vector<AcceptanceCondition::ptr>& acceptanceConditions) {
        // flattened vector storing unique acceptance sets
        auto numAcceptanceSets = 0;
        auto numStates = 0;
        for (const auto& acceptanceCondition : acceptanceConditions) {
            numAcceptanceSets += acceptanceCondition->getNumberOfAcceptanceSets();
        }
        uint64_t offset = 0;
        for (auto& ac : acceptanceConditions) {
            auto acceptanceExpr = cpphoafparser::addOffsetAcceptanceSets(ac->getAcceptanceExpression(), offset);
            auto liftedCondition = std::make_shared<AcceptanceCondition>(numStates, numAcceptanceSets, acceptanceExpr);

            for (uint64_t i = 0; i < ac->getNumberOfAcceptanceSets(); ++i) {
                liftedCondition->getAcceptanceSet(i+offset) = ac->getAcceptanceSet(i);
            }

            offset += ac->getNumberOfAcceptanceSets();
            ac = liftedCondition;
        }

        for (uint64_t i = 0; i < acceptanceConditions.size(); ++i) {
            for (uint64_t j = 0; j < acceptanceConditions.size(); ++j) {
                if (i == j) continue;
                for (uint64_t k = 0; k < acceptanceConditions[i]->getNumberOfAcceptanceSets(); ++k) {
                    if (!acceptanceConditions[i]->getAcceptanceSet(k).getNumberOfSetBits()) continue;
                    acceptanceConditions[j]->getAcceptanceSet(k) = acceptanceConditions[i]->getAcceptanceSet(k);
                }
            }
        }

        return acceptanceConditions;
    }

    static std::vector<AcceptanceCondition::ptr> getAllCombinations(std::vector<AcceptanceCondition::ptr> acceptanceConditions) {
        uint64_t numberAcceptanceConditionCombinations = pow(2, acceptanceConditions.size());
        std::vector<AcceptanceCondition::ptr> objectivesCombinations(numberAcceptanceConditionCombinations);

        acceptanceConditions = liftAcceptanceSets(acceptanceConditions);
        auto numAccSets = acceptanceConditions[0]->getNumberOfAcceptanceSets();
        auto numStates = acceptanceConditions[0]->getAcceptanceSet(0).size();

        for (uint64_t i = 0; i < numberAcceptanceConditionCombinations; i++) {
            auto acceptanceExpression = std::make_shared<cpphoafparser::BooleanExpression<cpphoafparser::AtomAcceptance>>(true);
            for (uint64_t j = 0; j < acceptanceConditions.size(); j++) {
                if (i & (1 << j)) {
                    acceptanceExpression = acceptanceConditions[j]->getAcceptanceExpression() & acceptanceExpression;
                } else {
                    acceptanceExpression = !acceptanceConditions[j]->getAcceptanceExpression() & acceptanceExpression;
                }
            }
            acceptanceExpression = toDNF(acceptanceExpression);
            STORM_LOG_INFO(acceptanceExpression->toString());
            objectivesCombinations[i] = std::make_shared<automata::AcceptanceCondition>(numStates, numAccSets, acceptanceExpression);

            // set acceptance sets
            for (uint j = 0; j < numAccSets; j++) {
                objectivesCombinations[i]->getAcceptanceSet(j) = acceptanceConditions[0]->getAcceptanceSet(j);
            }
        }

        return objectivesCombinations;
    }

private:
    std::vector<AcceptanceCondition::ptr> acceptanceConditions;

};
}
}
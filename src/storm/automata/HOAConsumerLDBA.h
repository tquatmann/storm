#pragma once

#include "storm/automata/APSet.h"
#include "storm/automata/HOAConsumerHeader.h"
#include "storm/exceptions/InvalidOperationException.h"
#include "storm/exceptions/OutOfRangeException.h"
#include "storm/solver/SmtSolver.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/expressions/ExpressionManager.h"
#include "storm/utility/solver.h"

#include "cpphoafparser/consumer/hoa_consumer.hh"
#include "cpphoafparser/util/implicit_edge_helper.hh"

#include <boost/optional.hpp>
#include <exception>

namespace storm {
namespace automata {

class HOAConsumerLDBA : public HOAConsumerHeader {
   private:
    AcceptanceCondition::ptr acceptance;

    // The consumer will build a non-deterministic automaton.
    LimitDeterministicAutomaton::ptr ldba;
    cpphoafparser::ImplicitEdgeHelper* helper = nullptr;

    std::shared_ptr<storm::expressions::ExpressionManager> expressionManager;
    std::vector<storm::expressions::Variable> apVariables;
    std::unique_ptr<storm::solver::SmtSolver> solver;

    storm::storage::BitVector seenEdges;

   public:
    typedef std::shared_ptr<HOAConsumerLDBA> ptr;

    HOAConsumerLDBA() : seenEdges(0) {
        expressionManager.reset(new storm::expressions::ExpressionManager());
        storm::utility::solver::SmtSolverFactory factory;
        solver = factory.create(*expressionManager);
    }

    ~HOAConsumerLDBA() {
        delete helper;
    }

    LimitDeterministicAutomaton::ptr getNDA() {
        return ldba;
    }

    void notifyBodyStart() override {
        if (!header.numberOfStates) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: Missing number-of-states header");
        }

        acceptance = header.getAcceptanceCondition();
        ldba.reset(new LimitDeterministicAutomaton(header.apSet, *header.numberOfStates, *header.startState, acceptance));

        helper = new cpphoafparser::ImplicitEdgeHelper(header.apSet.size());
        seenEdges.resize(*header.numberOfStates * helper->getEdgesPerState());

        for (const std::string& ap : header.apSet.getAPs()) {
            apVariables.push_back(expressionManager->declareBooleanVariable(ap));
        }
    }

    void addState(unsigned int id, std::shared_ptr<std::string> /*info*/, label_expr::ptr labelExpr, std::shared_ptr<int_list> accSignature) override {
        if (accSignature) {
            for (unsigned int accSet : *accSignature) {
                acceptance->getAcceptanceSet(accSet).set(id);
            }
        }

        if (labelExpr) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: State-labeled automata not supported");
        }

        helper->startOfState(id);
    }

    void addEdgeImplicit(unsigned int stateId, const int_list& conjSuccessors, std::shared_ptr<int_list> accSignature) override {
        std::size_t edgeIndex = helper->nextImplicitEdge();

        // Allow non-determinism by accepting multiple successors.
        if (conjSuccessors.empty()) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: Missing successor states in edge definition.");
        }

        if (accSignature) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: Does not support transition-based acceptance");
        }

        for (auto successor : conjSuccessors) {
            ldba->addSuccessor(stateId, edgeIndex, successor);
        }
        markEdgeAsSeen(stateId, edgeIndex);
    }

    void addEdgeWithLabel(unsigned int stateId, label_expr::ptr labelExpr, const int_list& conjSuccessors, std::shared_ptr<int_list> accSignature) override {
        // Allow non-determinism.
        if (conjSuccessors.empty()) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: Missing successor states in edge definition.");
        }

        if (accSignature) {
            throw std::runtime_error("Parsing non-deterministic HOA automaton: Does not support transition-based acceptance");
        }

        solver->reset();
        solver->add(labelToStormExpression(labelExpr));

        solver->allSat(apVariables, [this, stateId, &conjSuccessors](storm::expressions::SimpleValuation& valuation) {
            APSet::alphabet_element edgeIndex = header.apSet.elementAllFalse();
            for (std::size_t i = 0; i < apVariables.size(); i++) {
                if (valuation.getBooleanValue(apVariables[i])) {
                    edgeIndex = header.apSet.elementAddAP(edgeIndex, i);
                }
            }

            // Add all successors for the current label.
            for (auto successor : conjSuccessors) {
                ldba->addSuccessor(stateId, edgeIndex, successor);
            }
            markEdgeAsSeen(stateId, edgeIndex);

            return true;
        });
    }

    void notifyEndOfState(unsigned int /*stateId*/) override {
        helper->endOfState();
    }

    void notifyEnd() override {
        STORM_LOG_THROW(seenEdges.full(), storm::exceptions::InvalidOperationException, "HOA automaton has mismatch in number of edges, not complete?");
    }

    void notifyAbort() override {
        throw std::runtime_error("Parsing non-deterministic automaton: Automaton is incomplete (abort)");
    }

    void notifyWarning(const std::string& /*warning*/) override {
        // IGNORE
    }

   private:
    storm::expressions::Expression labelToStormExpression(label_expr::ptr labelExpr) {
        switch (labelExpr->getType()) {
            case label_expr::EXP_AND:
                return labelToStormExpression(labelExpr->getLeft()) && labelToStormExpression(labelExpr->getRight());
            case label_expr::EXP_OR:
                return labelToStormExpression(labelExpr->getLeft()) || labelToStormExpression(labelExpr->getRight());
            case label_expr::EXP_NOT:
                return !labelToStormExpression(labelExpr->getLeft());
            case label_expr::EXP_TRUE:
                return expressionManager->boolean(true);
            case label_expr::EXP_FALSE:
                return expressionManager->boolean(false);
            case label_expr::EXP_ATOM: {
                unsigned int apIndex = labelExpr->getAtom().getAPIndex();
                STORM_LOG_THROW(apIndex < apVariables.size(), storm::exceptions::OutOfRangeException,
                                "HOA automaton refers to non-existing atomic proposition");
                return apVariables.at(apIndex).getExpression();
            }
        }
        throw std::runtime_error("Unknown label expression operator");
    }

    void markEdgeAsSeen(std::size_t stateId, std::size_t edgeIndex) {
        seenEdges.set(stateId * helper->getEdgesPerState() + edgeIndex);
    }
};

}  // namespace automata
}  // namespace storm
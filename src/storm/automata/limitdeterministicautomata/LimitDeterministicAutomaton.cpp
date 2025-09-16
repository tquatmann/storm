#include "LimitDeterministicAutomaton.h"

#include <storm/automata/APSet.h>

#include "HOAConsumerLDBA.h"
#include "cpphoafparser/consumer/hoa_intermediate_check_validity.hh"
#include "cpphoafparser/parser/hoa_parser.hh"
#include "cpphoafparser/parser/hoa_parser_helper.hh"

#include "storm/automata/AcceptanceCondition.h"
#include "storm/exceptions/FileIoException.h"
#include "storm/io/file.h"
#include "storm/utility/macros.h"

namespace storm {
namespace automata {

LimitDeterministicAutomaton::LimitDeterministicAutomaton(APSet apSet, storm::storage::sparse::state_type numberOfStates,
                                                         storm::storage::sparse::state_type initialState, std::shared_ptr<AcceptanceCondition> acceptance)
    : apSet(apSet), numberOfStates(numberOfStates), initialState(initialState), acceptance(acceptance) {
    edgesPerState = apSet.alphabetSize();
    successors.resize(numberOfStates * edgesPerState, storage::BitVector(numberOfStates));
    predecessors.resize(numberOfStates, storage::BitVector(numberOfStates));
}

const APSet& LimitDeterministicAutomaton::getAPSet() const {
    return apSet;
}

storm::storage::sparse::state_type LimitDeterministicAutomaton::getInitialState() const {
    return initialState;
}

const storm::storage::BitVector& LimitDeterministicAutomaton::getSuccessors(storm::storage::sparse::state_type from, APSet::alphabet_element label) const {
    storm::storage::sparse::state_type index = from * edgesPerState + label;
    return successors.at(index);
}

void LimitDeterministicAutomaton::addSuccessor(storm::storage::sparse::state_type from, APSet::alphabet_element label,
                                               storm::storage::sparse::state_type successor) {
    storm::storage::sparse::state_type index = from * edgesPerState + label;
    successors.at(index).set(successor);
    predecessors.at(successor).set(from);
}

storm::storage::sparse::state_type LimitDeterministicAutomaton::getNumberOfStates() const {
    return numberOfStates;
}

storm::storage::sparse::state_type LimitDeterministicAutomaton::getNumberOfEdgesPerState() const {
    return edgesPerState;
}

std::shared_ptr<AcceptanceCondition> LimitDeterministicAutomaton::getAcceptance() const {
    return acceptance;
}

storm::storage::BitVector LimitDeterministicAutomaton::getAcceptingPart() const {
    auto expr = acceptance->getAcceptanceExpression();
    STORM_LOG_ASSERT(expr->isAtom() && expr->getAtom().getType() == cpphoafparser::AtomAcceptance::TEMPORAL_INF && acceptance->getNumberOfAcceptanceSets() == 1,
                     "Büchi acceptance condition has to be of the form INF(0)");

    auto acceptingPart = acceptance->getAcceptanceSet(0);

    std::queue<uint64_t> todo;
    for (auto const& state : acceptingPart) {
        todo.push(state);
    }

    while (!todo.empty()) {
        auto state = todo.front();
        todo.pop();

        for (auto const& next : predecessors.at(state) & ~acceptingPart) {
            todo.push(next);
        }
        acceptingPart |= predecessors.at(state);
    }

    return acceptingPart;
}

void LimitDeterministicAutomaton::printHOA(std::ostream& out) const {
    out << "HOA: v1\n";

    out << "States: " << numberOfStates << "\n";

    out << "Start: " << initialState << "\n";

    out << "AP: " << apSet.size();
    for (unsigned int i = 0; i < apSet.size(); i++) {
        out << " " << cpphoafparser::HOAParserHelper::quote(apSet.getAP(i));
    }
    out << "\n";

    out << "Acceptance: " << acceptance->getNumberOfAcceptanceSets() << " " << *acceptance->getAcceptanceExpression() << "\n";

    out << "--BODY--" << "\n";

    for (storm::storage::sparse::state_type s = 0; s < getNumberOfStates(); s++) {
        out << "State: " << s;
        out << " {";
        bool first = true;
        for (unsigned int i = 0; i < acceptance->getNumberOfAcceptanceSets(); i++) {
            if (acceptance->getAcceptanceSet(i).get(s)) {
                if (!first)
                    out << " ";
                first = false;
                out << i;
            }
        }
        out << "}\n";
        for (storm::storage::sparse::state_type label = 0; label < getNumberOfEdgesPerState(); label++) {
            const auto& succs = getSuccessors(s, label);
            if (succs.empty()) {
                // This case should not happen in a valid automaton.
                // Every state must have a successor for every label.
            } else {
                out << "(";
                for (size_t i = 0; i < succs.size(); ++i) {
                    out << succs[i] << (i == succs.size() - 1 ? "" : " & ");
                }
                out << ")\n";
            }
        }
    }
}

LimitDeterministicAutomaton::ptr LimitDeterministicAutomaton::parse(std::istream& in) {
    HOAConsumerLDBA::ptr consumer(new HOAConsumerLDBA());
    cpphoafparser::HOAIntermediateCheckValidity::ptr validator(new cpphoafparser::HOAIntermediateCheckValidity(consumer));
    cpphoafparser::HOAParser::parse(in, validator);

    return consumer->getNDA();
}

LimitDeterministicAutomaton::ptr LimitDeterministicAutomaton::parseFromFile(const std::string& filename) {
    std::ifstream in;
    storm::io::openFile(filename, in);
    auto ldba = parse(in);
    storm::io::closeFile(in);

    STORM_LOG_INFO("Nondeterministic automaton from HOA file '" << filename << "' has " << ldba->getNumberOfStates() << " states, " << ldba->getAPSet().size()
                                                                << " atomic propositions and " << *ldba->getAcceptance()->getAcceptanceExpression()
                                                                << " as acceptance condition.");
    return ldba;
}

}  // namespace automata
}  // namespace storm
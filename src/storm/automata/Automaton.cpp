#include "storm/automata/Automaton.h"

#include "cpphoafparser/consumer/hoa_intermediate_check_validity.hh"
#include "cpphoafparser/parser/hoa_parser.hh"
#include "cpphoafparser/parser/hoa_parser_helper.hh"

#include "storm/automata/AcceptanceCondition.h"
#include "storm/automata/HOAConsumerDA.h"
#include "storm/automata/HOAConsumerLDBA.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/io/file.h"
#include "storm/utility/macros.h"

namespace storm::automata {

template<AutomatonType Type>
Automaton<Type>::Automaton(APSet apSet, std::size_t numberOfStates, StateIndex initialState, AcceptanceCondition::ptr acceptance)
    : apSet(std::move(apSet)), numberOfStates(numberOfStates), initialState(initialState), acceptance(acceptance) {
    edgesPerState = this->apSet.alphabetSize();
    SuccessorType const defaultSuccessor{numberOfStates};  // Either an invalid index or a BitVector of the right size.
    successors.resize(numberOfStates * edgesPerState, defaultSuccessor);
}

template<AutomatonType Type>
typename Automaton<Type>::StateIndex Automaton<Type>::getInitialState() const {
    return initialState;
}

template<AutomatonType Type>
const APSet& Automaton<Type>::getAPSet() const {
    return apSet;
}

template<AutomatonType Type>
typename Automaton<Type>::StateIndex Automaton<Type>::getSuccessor(StateIndex from, APSet::alphabet_element label) const
    requires isDeterministic
{
    auto const index = from * edgesPerState + label;
    STORM_LOG_ASSERT(successors.at(index) < numberOfStates, "The automaton is not complete, missing successor for state " << from << " and label " << label);
    return successors.at(index);
}

template<AutomatonType Type>
void Automaton<Type>::setSuccessor(StateIndex from, APSet::alphabet_element label, StateIndex successor)
    requires isDeterministic
{
    STORM_LOG_ASSERT(successor < numberOfStates, "The successor " << successor << " is not a valid state index.");
    auto const index = from * edgesPerState + label;
    successors.at(index) = successor;
}

template<AutomatonType Type>
typename Automaton<Type>::SuccessorType Automaton<Type>::getSuccessors(StateIndex from, APSet::alphabet_element label) const
    requires(!isDeterministic)
{
    auto const index = from * edgesPerState + label;
    return successors.at(index);
}

template<AutomatonType Type>
void Automaton<Type>::addSuccessor(StateIndex from, APSet::alphabet_element label, StateIndex successor) {
    STORM_LOG_ASSERT(successor < numberOfStates, "The successor " << successor << " is not a valid state index.");
    auto const index = from * edgesPerState + label;
    auto& succ = successors.at(index);
    if constexpr (isDeterministic) {
        STORM_LOG_ASSERT(succ == numberOfStates, "The automaton is not deterministic, multiple successors for state " << from << " and label " << label);
        succ = successor;
    } else {
        succ.set(successor);
    }
}

template<AutomatonType Type>
std::size_t Automaton<Type>::getNumberOfStates() const {
    return numberOfStates;
}

template<AutomatonType Type>
std::size_t Automaton<Type>::getNumberOfEdgesPerState() const {
    return edgesPerState;
}

template<AutomatonType Type>
AcceptanceCondition::ptr Automaton<Type>::getAcceptance() const {
    return acceptance;
}

template<AutomatonType Type>
storm::storage::BitVector Automaton<Type>::getAcceptingPart() const
    requires(Type == AutomatonType::LDBA)
{
    auto expr = acceptance->getAcceptanceExpression();
    STORM_LOG_ASSERT(expr->isAtom() && expr->getAtom().getType() == cpphoafparser::AtomAcceptance::TEMPORAL_INF && acceptance->getNumberOfAcceptanceSets() == 1,
                     "Büchi acceptance condition has to be of the form INF(0)");
    auto acceptingPart = acceptance->getAcceptanceSet(0);

    // Find all states reachable from an accepting state
    std::vector<StateIndex> queue(acceptingPart.begin(), acceptingPart.end());
    while (!queue.empty()) {
        auto const state = queue.back();
        queue.pop_back();

        for (uint64_t label = 0; label < edgesPerState; label++) {
            STORM_LOG_ASSERT(getSuccessors(state, label).hasUniqueSetBit(), "The automaton is not limit deterministic: nondeterministic choice for state "
                                                                                << state << " and label " << label << " reachable from an accepting state.");
            for (auto const next : getSuccessors(state, label)) {
                if (!acceptingPart.get(next)) {
                    queue.push_back(next);
                    acceptingPart.set(next);
                }
            }
        }
    }
    return acceptingPart;
}

template<AutomatonType Type>
void Automaton<Type>::printHOA(std::ostream& out) const {
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

    for (StateIndex s = 0; s < getNumberOfStates(); s++) {
        out << "State: " << s;
        out << " {";
        bool first = true;
        for (uint64_t i = 0; i < acceptance->getNumberOfAcceptanceSets(); i++) {
            if (acceptance->getAcceptanceSet(i).get(s)) {
                if (!first)
                    out << " ";
                first = false;
                out << i;
            }
        }
        out << "}\n";
        for (uint64_t label = 0; label < edgesPerState; label++) {
            if constexpr (isDeterministic) {
                out << getSuccessor(s, label) << "\n";
            } else {
                const auto& succs = getSuccessors(s, label);
                STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "HOA output of nondeterministic automata is not yet implemented");
                STORM_LOG_ASSERT(!succs.empty(), "The automaton is not complete, missing successor for state " << s << " and label " << label);
            }
        }
    }
}

template<AutomatonType Type>
typename Automaton<Type>::ptr Automaton<Type>::parse(std::istream& in) {
    using ConsumerType = std::conditional_t<isDeterministic, HOAConsumerDA, HOAConsumerLDBA>;
    typename ConsumerType::ptr consumer(new ConsumerType());
    cpphoafparser::HOAIntermediateCheckValidity::ptr validator(new cpphoafparser::HOAIntermediateCheckValidity(consumer));
    cpphoafparser::HOAParser::parse(in, validator);

    if constexpr (isDeterministic) {
        return consumer->getDA();
    } else {
        return consumer->getNDA();
    }
}

template<AutomatonType Type>
Automaton<Type>::ptr Automaton<Type>::parseFromFile(const std::string& filename) {
    std::ifstream in;
    storm::io::openFile(filename, in);
    auto aut = parse(in);
    storm::io::closeFile(in);

    STORM_LOG_INFO(Type << " Automaton from HOA file '" << filename << "' has " << aut->getNumberOfStates() << " states, " << aut->getAPSet().size()
                        << " atomic propositions and " << *aut->getAcceptance()->getAcceptanceExpression() << " as acceptance condition.");
    return aut;
}

template class Automaton<AutomatonType::DA>;
template class Automaton<AutomatonType::LDBA>;

}  // namespace storm::automata

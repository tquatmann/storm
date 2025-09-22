#pragma once

#include "storm/automata/AutomatonType.h"

namespace storm::automata {

template<AutomatonType Type>
class Automaton;

using DeterministicAutomaton = Automaton<AutomatonType::DA>;
using LimitDeterministicAutomaton = Automaton<AutomatonType::LDBA>;

}  // namespace storm::automata

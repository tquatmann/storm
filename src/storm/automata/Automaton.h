#pragma once

#include <iostream>
#include <memory>
#include <span>

#include "storm/automata/APSet.h"
#include "storm/automata/AutomatonForward.h"

namespace storm {

namespace storage {
class BitVector;
}

namespace automata {
class AcceptanceCondition;

template<AutomatonType Type>
class Automaton {
   public:
    using ptr = std::shared_ptr<Automaton<Type>>;
    using StateIndex = uint64_t;
    static constexpr AutomatonType type = Type;
    static constexpr bool isDeterministic = (Type == AutomatonType::DA);
    using SuccessorType = std::conditional_t<isDeterministic, StateIndex, storm::storage::BitVector>;

    Automaton(APSet apSet, std::size_t numberOfStates, StateIndex initialState, std::shared_ptr<AcceptanceCondition> acceptance);

    APSet const& getAPSet() const;
    StateIndex getInitialState() const;

    // Successors
    StateIndex getSuccessor(StateIndex from, APSet::alphabet_element label) const
        requires isDeterministic;
    void setSuccessor(StateIndex from, APSet::alphabet_element label, StateIndex successor)
        requires isDeterministic;
    SuccessorType getSuccessors(StateIndex from, APSet::alphabet_element label) const
        requires(!isDeterministic);
    void addSuccessor(StateIndex from, APSet::alphabet_element label, StateIndex successor);

    std::size_t getNumberOfStates() const;
    std::size_t getNumberOfEdgesPerState() const;

    std::shared_ptr<AcceptanceCondition> getAcceptance() const;

    /*!
     * Computes the partition of the LDBA states into initial and accepting part.
     * The accepting part contains all accepting states and all of the contained states are deterministic.
     * Furthermore, the accepting part can only be entered via a nondeterministic state.
     * @return the set of automaton states that are in the accepting part.
     */
    storm::storage::BitVector computeAcceptingPart() const
        requires(Type == AutomatonType::LDBA);

    void printHOA(std::ostream& out) const;

    static ptr parse(std::istream& in);
    static ptr parseFromFile(const std::string& filename);

   private:
    APSet apSet;
    std::size_t numberOfStates;
    StateIndex initialState;
    std::size_t edgesPerState;
    std::shared_ptr<AcceptanceCondition> acceptance;
    std::vector<SuccessorType> successors;
};

}  // namespace automata
}  // namespace storm

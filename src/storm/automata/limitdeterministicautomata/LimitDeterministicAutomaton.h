#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include "storm/automata/APSet.h"
#include "storm/automata/AcceptanceCondition.h"

namespace storm {
namespace automata {

class LimitDeterministicAutomaton {
   public:
    using ptr = std::shared_ptr<LimitDeterministicAutomaton>;

    LimitDeterministicAutomaton(APSet apSet, storm::storage::sparse::state_type numberOfStates, storm::storage::sparse::state_type initialState,
                                std::shared_ptr<AcceptanceCondition> acceptance);

    const APSet& getAPSet() const;
    storm::storage::sparse::state_type getInitialState() const;
    const storm::storage::BitVector& getSuccessors(storm::storage::sparse::state_type from, APSet::alphabet_element label) const;
    void addSuccessor(storm::storage::sparse::state_type from, APSet::alphabet_element label, storm::storage::sparse::state_type successor);

    storm::storage::sparse::state_type getNumberOfStates() const;
    storm::storage::sparse::state_type getNumberOfEdgesPerState() const;

    std::shared_ptr<AcceptanceCondition> getAcceptance() const;
    storm::storage::BitVector getAcceptingPart() const;

    void printHOA(std::ostream& out) const;

    static LimitDeterministicAutomaton::ptr parse(std::istream& in);
    static LimitDeterministicAutomaton::ptr parseFromFile(const std::string& filename);

   private:
    APSet apSet;
    storm::storage::sparse::state_type numberOfStates;
    storm::storage::sparse::state_type initialState;
    storm::storage::sparse::state_type edgesPerState;
    std::shared_ptr<AcceptanceCondition> acceptance;
    std::vector<storm::storage::BitVector> successors;
    std::vector<storm::storage::BitVector> predecessors;
};

}  // namespace automata
}  // namespace storm
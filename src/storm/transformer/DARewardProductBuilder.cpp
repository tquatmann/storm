#include "DARewardProductBuilder.h"
#include "../storage/MaximalEndComponentDecomposition.h"

//#define DA_READY

namespace storm {
namespace transformer {
template<typename ValueType>
DAProduct<Model>::ptr DARewardProductBuilder<Model>::build(const storm::storage::SparseMatrix<typename Model::ValueType>& originalMatrix, const storm::storage::BitVector& statesOfInterest) {
    transformer::DAProductBuilder productBuilder(da, statesForAP);
    auto product = productBuilder.build<Model>(originalMatrix, statesOfInterest);

    /*
    STORM_LOG_INFO("Computing accepting states for acceptance condition " << *acceptance.getAcceptanceExpression());
    if (acceptance.getAcceptanceExpression()->isTRUE()) {
        STORM_LOG_INFO(" TRUE -> all states accepting (assumes no deadlock in the model)");
        return storm::storage::BitVector(transitionMatrix.getRowGroupCount(), true);
    } else if (acceptance.getAcceptanceExpression()->isFALSE()) {
        STORM_LOG_INFO(" FALSE -> all states rejecting");
        return storm::storage::BitVector(transitionMatrix.getRowGroupCount(), false);
    }
    */

    return product;
}

#ifdef DA_READY

template<typename Model>
void DARewardProductBuilder<Model>::modify(storm::storage::SparseMatrix<typename Model::ValueType> const& transitionMatrix,
                                            storm::storage::SparseMatrix<typename Model::ValueType> const& backwardTransitions) {
    // MatrixBuilder to build the transition matrix for the demerged MDP
    storage::SparseMatrixBuilder<typename Model::ValueType> transitionMatrixBuilder(0, 0, 0, false, true, false, 0);

    //Compute MECs for entire MDP
    storm::storage::MaximalEndComponentDecomposition<Model> mecs(transitionMatrix, backwardTransitions);
    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), std::numeric_limits<uint64_t>::max());
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        for (auto const& [state, _] : mec) {
            stateToMec[state] = mec_counter;
        }
        ++mec_counter;
    }
    uint64_t InvalidMECIndex = std::numeric_limits<uint64_t>::max();

    //rows to be added, for states representing the MECs
    storage::SparseMatrixBuilder<typename Model::ValueType> representativeStateMatrixBuilder(0, 0, 0, false, true, false, 0);
    //im zweiten schritt einfach klopierend und mac actions zwischendurch hinzufuegen


    for (auto const& state: transitionMatrix.getRowGroupCount()) {
        uint64_t const mecIndex = stateToMec[state];
        transitionMatrix.newRowGroup();

        for (auto const& choice: transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);

            if (mecIndex == InvalidMECIndex) {
                modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size(), transitionMatrix.getRowGroupCount());
            }

            bool choiceIsInMec = mecs[mecIndex].getChoicesForState(state).find(choice) != mecs[mecIndex].getChoicesForState(state).end();

            if (mecIndex != InvalidMECIndex && choiceIsInMec) {
                // copy values from
                for (auto& entry: row) transitionMatrixBuilder.addNextValue(choice, entry.getColumn(), entry.getValue());
            }

            if (mecIndex == InvalidMECIndex && !choiceIsInMec) {
                // actions are unique
                modifyStateActionPair(representativeStateMatrixBuilder, stateToMec, row, mecs.size(), transitionMatrix.getRowGroupCount());
            }
        }
    }

    auto representativeStateMatrix = representativeStateMatrixBuilder.build();

    // compute MACs
    // returns list of macs
    std::list<storage::MaximalEndComponent> accEcs = computeAcceptingEcs(acceptance, transitionMatrix, backwardTransitions);

    // add MEC matrix to main matrix
    std::vector<std::list<uint64_t> > mecToMacStates(mecs.size(), std::list<uint64_t>());

    for (auto const& mac: accEcs) {
        uint64_t representativeEcState = mac.begin()->first;
        auto const mec_index = stateToMec[representativeEcState];
        STORM_LOG_ASSERT(mec_index < mecs.size(), "MEC index out of range.");
        mecToMacStates[mec_index].push_back(representativeEcState);
    }

    for (auto const& state: representativeStateMatrix.getRowGroupIndices()) {
        transitionMatrix.newRowGroup();
        for (auto const& choice: representativeStateMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);
            for (auto const& entry: row) {
                transitionMatrixBuilder.addNextValue(numRows, entry.getColumn(), entry.getValue());
            }
        }

        for (auto const& macState: mecToMacStates[state]) {
            transitionMatrixBuilder.addNextValue(numRows, macState + state_offset, 1);
        }
    }

    //STORM_LOG_INFO("Found " << accMECs << " accepting MECs (considered " << allMECs << " MECs).");


    // add MACs to matrix
    for (auto const& mac: accEcs) {
        //states and choices from original matrix
        for (auto const& [state, choices]: mac) {
            transitionMatrixBuilder.newRowGroup();

            for (auto const& choice: choices) {
                auto row = transitionMatrix.getRow(choice);
                for (auto& entry: row) transitionMatrixBuilder.addNextValue(choice + choice_offset, entry.getColumn() + state_offset, entry.getValue());
            }
        }
    }

    return transitionMatrixBuilder.build();
}

template<typename Model>
void DARewardProductBuilder<Model>::modifyStateActionPair(storage::SparseMatrixBuilder<typename Model::ValueType>& builder,
                                                                                            std::vector<uint64_t>& stateToMec,
                                                                                            row,
                                                                                            uint64_t numMecs,
                                                                                            uint64_t numStates) {
    std::vector<Model::ValueType> mecValues(numMecs, 0);
    uint64_t choice = builder.getLastRow() + 1;
    uint64_t InvalidMECIndex = std::numeric_limits<uint64_t>::max();

    for (auto& entry: row) {
        uint64_t mecIndexNextState = stateToMec[entry.getColumn()];

        if (mecIndexNextState != InvalidMECIndex) {
            builder.addNextValue(choice, entry.getColumn(), entry.getValue());
        }
        else {
            mecValues[mecIndexNextState] += entry.getValue();
        }
    }

    for (uint64_t i=0; auto& entry: mecValues) {
        if (entry) {
            builder.addNextValue(choice, numStates + i,entry);
        }
        i++;
    }
}

template<typename ValueType>
std::list<storage::MaximalEndComponent> DARewardProductBuilder<ValueType>::computeAcceptingECs(automata::AcceptanceCondition const& acceptance,
                                                                                            storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                                            storm::storage::SparseMatrix<ValueType> const& backwardTransitions) {
    // list of maximum accepting end components in the product MDP
    std::list<storage::MaximalEndComponent> macs;
    std::vector<std::vector<automata::AcceptanceCondition::acceptance_expr::ptr>> dnf = acceptance.extractFromDNF();

    for (auto const& conjunction : dnf) {
        // get the states of the mdp that are on a MEC and don't violate Fins of the conjunction
        storm::storage::BitVector allowed = transitionMatrix.getRowGroupCount(); //maybe late start with mecStates

        for (auto const& literal : conjunction) {
            if (literal->isTRUE()) {
                // skip
            } else if (literal->isFALSE()) {
                allowed.clear();
                break;
            } else if (literal->isAtom()) {
                const cpphoafparser::AtomAcceptance& atom = literal->getAtom();
                if (atom.getType() == cpphoafparser::AtomAcceptance::TEMPORAL_FIN) {
                    // only deal with FIN, ignore INF here
                    const storm::storage::BitVector& accSet = acceptance.getAcceptanceSet(atom.getAcceptanceSet());
                    if (atom.isNegated()) {
                        // allowed = allowed \ ~accSet = allowed & accSet
                        allowed &= accSet;
                    } else {
                        // allowed = allowed \ accSet = allowed & ~accSet
                        allowed &= ~accSet;
                    }
                }
            }
        }

        if (allowed.empty()) {
            // skip
            continue;
        }

        // Compute MECs in the allowed fragment
        storm::storage::MaximalEndComponentDecomposition<ValueType> allowedECs(transitionMatrix, backwardTransitions, allowed);
        for (const auto& ec : allowedECs) {
            bool accepting = true;
            for (auto const& literal : conjunction) {
                if (literal->isTRUE()) {
                    // skip

                } else if (literal->isFALSE()) {
                    accepting = false;
                    break;
                } else if (literal->isAtom()) {
                    const cpphoafparser::AtomAcceptance& atom = literal->getAtom();
                    const storm::storage::BitVector& accSet = acceptance.getAcceptanceSet(atom.getAcceptanceSet());
                    if (atom.getType() == cpphoafparser::AtomAcceptance::TEMPORAL_INF) {
                        if (atom.isNegated()) {
                            if (!ec.containsAnyState(~accSet)) {
                                accepting = false;
                                break;
                            }
                        } else {
                            if (!ec.containsAnyState(accSet)) {
                                accepting = false;
                                break;
                            }
                        }

                    } else if (atom.getType() == cpphoafparser::AtomAcceptance::TEMPORAL_FIN) {
                        // Do only sanity checks here.
                        STORM_LOG_ASSERT(atom.isNegated() ? !ec.containsAnyState(~accSet) : !ec.containsAnyState(accSet),
                                         "EC contains Fin-states, which should have been removed");
                    }
                }
            }

            if (accepting) {
                macs.push_back(ec);
            }
        }
    }

    return macs;
}

#endif
}
}


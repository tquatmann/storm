#include "DARewardProductBuilder.h"

#include "../modelchecker/CheckTask.h"
#include "../modelchecker/helper/ltl/SparseLTLHelper.h"
#include "../models/sparse/MarkovAutomaton.h"
#include "../storage/prism/RewardModel.h"
#include "storm/automata/DeterministicAutomaton.h"
#include "storm/automata/LTL2DeterministicAutomaton.h"
#include "storm/logic/ExtractMaximalStateFormulasVisitor.h"

#include <gmpxx.h>  //???
#include "storm/exceptions/InvalidPropertyException.h"

#define DEBUG_VAR(x) std::cout << #x << " = " << x << std::endl;

namespace storm {
namespace transformer {

template <typename T>
void printMatrix(storm::storage::SparseMatrix<T> matrix);

template<typename ValueType, typename RewardModelType>
std::shared_ptr<DARewardProduct<ValueType>> DARewardProductBuilder<ValueType, RewardModelType>::build() {
    auto transitionMatrix = product.getProductModel().getTransitionMatrix();
    auto backwardTransitions = product.getProductModel().getBackwardTransitions();
    auto acceptance = *product.getAcceptance();

    STORM_LOG_INFO("Building transition matrix...");

    // map every state/action in the modified model to the original state/action in the product model
    std::vector<uint64_t> stateToProductState, choiceToModelChoice;
    std::list<uint64_t> reachingAccEcChoices;

    // MatrixBuilder to build the transition matrix for the demerged MDP
    storage::SparseMatrixBuilder<ValueType> transitionMatrixBuilder(0, 0, 0, false, true, 0);
    auto lastRow = [&transitionMatrixBuilder](){ return transitionMatrixBuilder.getLastRow();};

    //Compute MECs for entire MDP
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(transitionMatrix, backwardTransitions);
    uint64_t totalNumberStates = transitionMatrix.getRowGroupCount() + mecs.size();
    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), InvalidIndex);
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        std::cout << "MEC contains:" << std::endl;
        for (auto const& [state, _] : mec) {
            stateToMec[state] = mec_counter;
        }
        ++mec_counter;
    }

    /*
    // set initial states of modified model
    storm::storage::BitVector initialStatesModified(totalNumberStates, false);
    for (uint64_t prodInitState: initialStates) {
        uint64_t mecIndex = stateToMec[prodInitState];
        if (mecIndex == InvalidIndex) {
            initialStatesModified.set(prodInitState, true);
        } else {
            initialStatesModified.set(transitionMatrix.getRowGroupCount() + mecIndex, true);
        }
    }
    */

    // compute the accepting end components
    std::list<storage::MaximalEndComponent> accEcs = computeAcceptingECs(acceptance, transitionMatrix, backwardTransitions);
    for (auto& mac: accEcs) {
        totalNumberStates += mac.size();
    }
    stateToProductState.resize(totalNumberStates);

    // A mapping of each mec to the actions that are enabled in a state in it but leave the mec with some probability
    std::vector<std::list<uint64_t>> mecsToLeavingActions(mecs.size(), std::list<uint64_t>());
    bool isMatrixBuilderEmpty = true;

    // Iterates over all state-action pairs and modifies them for new model
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        uint64_t const mecIndex = stateToMec[state];
        uint64_t startRowGroup = isMatrixBuilderEmpty ? 0 : lastRow() + 1;
        transitionMatrixBuilder.newRowGroup(startRowGroup);
        stateToProductState[state] = state;

        for (auto const& choice: transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);
            uint64_t nextRow = isMatrixBuilderEmpty ? 0 : lastRow() + 1;

            if (mecIndex == InvalidIndex) {
                modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size(), isMatrixBuilderEmpty);
            }

            bool choiceIsInMec = (mecIndex != InvalidIndex) && mecs[mecIndex].containsChoice(state, choice);

            // The choice is part of an MEC and therefore cannot leave it
            if (choiceIsInMec) {
                for (auto& entry: row) {
                    transitionMatrixBuilder.addNextValue(nextRow, entry.getColumn(), entry.getValue());
                }
            }

            if (mecIndex != InvalidIndex && !choiceIsInMec) {
                // actions are unique
                mecsToLeavingActions[mecIndex].push_back(choice);
            }

            isMatrixBuilderEmpty = false;
        }
    }
    std::cout << "MECs: " << mecsToLeavingActions.size() << std::endl;

    // add MEC matrix to main matrix
    std::vector<std::list<uint64_t> > mecToMacStates(mecs.size(), std::list<uint64_t>());
    // represents the first state of each MAC so we can later add an action that leads there a.s.
    uint64_t representativeColumn = transitionMatrix.getColumnCount() + mecs.size();

    // determining to which MEC the end component belongs to
    for (auto const& mac: accEcs) {
        uint64_t representativeEcState = mac.begin()->first;
        auto const mec_index = stateToMec[representativeEcState];
        STORM_LOG_ASSERT(mec_index < mecs.size(), "MEC index out of range.");
        mecToMacStates[mec_index].push_back(representativeColumn);
        representativeColumn += accEcs.size();
    }

    // add states representing MECs and the enabled actions that leave them
    for (uint64_t mecIndex = 0; mecIndex < mecs.size(); ++mecIndex) {
        transitionMatrixBuilder.newRowGroup(lastRow() + 1);
        // these states have no counterpart in the product model
        stateToProductState[mecIndex + transitionMatrix.getRowGroupCount()] = InvalidIndex;

        for (auto const& choice: mecsToLeavingActions[mecIndex]) {
            auto row = transitionMatrix.getRow(choice);
            modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size());
        }

        // actions that reach copies of accepting end components
        for (auto const& macState: mecToMacStates[mecIndex]) {
            uint64_t choice = lastRow() + 1;
            reachingAccEcChoices.push_back(choice);
            transitionMatrixBuilder.addNextValue(choice, macState, 1);
        }

        // action to representative state of MEC
        auto representativeState = mecs[mecIndex].begin()->first;
        transitionMatrixBuilder.addNextValue(lastRow() + 1, representativeState, 1);
    }

    uint64_t column_offset = transitionMatrix.getColumnCount() + mecs.size();
    uint64_t stateCounter = transitionMatrix.getRowGroupCount() + mecs.size();
    // add MACs to matrix
    for (auto const& mac: accEcs) {
        auto stateSet = mac.getStateSet();

        //states and choices from original matrix
        for (auto const& [state, choices]: mac) {
            transitionMatrixBuilder.newRowGroup(lastRow() + 1);
            stateToProductState[stateCounter++] = state;

            for (auto const& choice: choices) {
                auto row = transitionMatrix.getRow(choice);
                uint64_t choiceModified = lastRow() + 1;

                for (auto& entry: row) {
                    uint64_t nextState = std::distance(stateSet.begin(), stateSet.find(entry.getColumn()));
                    transitionMatrixBuilder.addNextValue(choiceModified, nextState + column_offset, entry.getValue());
                }
            }
        }
        column_offset += stateSet.size();
    }

    /*
    // replace product state by model state in mapping
    std::vector<uint64_t> stateToModelState(result.stateToModelState.size());
    for (uint64_t state = 0; state < result.stateToModelState.size(); state++) {
        uint64_t prodState = result.stateToModelState[state];

        if (prodState != InvalidIndex) {
            stateToModelState[state] = product->getModelState(prodState);
        } else {
            stateToModelState[state] = InvalidIndex;
        }
    }
    */

    return nullptr; // std::make_shared<DARewardProduct<ValueType>>(transitionMatrix);
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, std::vector<uint64_t>& stateToMec, Row row, uint64_t numMecs, bool builderEmpty) {
    std::vector<ValueType> mecValues(numMecs, 0);
    uint64_t choice = builderEmpty ? 0 : builder.getLastRow() + 1;
    ValueType sumEntries = 0;

    for (auto& entry: row) {
        uint64_t mecIndexNextState = stateToMec[entry.getColumn()];

        if (mecIndexNextState == std::numeric_limits<uint64_t>::max()) {
            builder.addNextValue(choice, entry.getColumn(), entry.getValue());
            sumEntries += entry.getValue();
        }
        else {
            mecValues[mecIndexNextState] += entry.getValue();
        }
    }

    for (uint64_t i=0; auto& entry: mecValues) {
        if (entry) {
            builder.addNextValue(choice, stateToMec.size() + i,entry);
            sumEntries += entry;
        }
        i++;
    }
}

template<typename ValueType, typename RewardModelType>
std::list<storage::MaximalEndComponent> DARewardProductBuilder<ValueType, RewardModelType>::computeAcceptingECs(automata::AcceptanceCondition const& acceptance,
                                                                                            storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                                            storm::storage::SparseMatrix<ValueType> const& backwardTransitions) {
    // list of maximum accepting end components in the product MDP
    std::list<storage::MaximalEndComponent> macs;
    std::vector<std::vector<automata::AcceptanceCondition::acceptance_expr::ptr>> dnf = acceptance.extractFromDNF();

    for (auto const& conjunction : dnf) {
        // get the states of the mdp that are on a MEC and don't violate Fins of the conjunction
        storm::storage::BitVector allowed(transitionMatrix.getRowGroupCount(), true); //maybe late start with mecStates

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

template <typename T>
void printMatrix(storage::SparseMatrix<T> matrix) {
    std::cout << "Rows: " << matrix.getRowCount() << ", Columns: " << matrix.getColumnCount() << ", Num entries: " << matrix.getEntryCount() << std::endl;

    for (uint64_t state = 0; state < matrix.getRowGroupCount(); state++) {
        std::cout << state << ": " << std::endl;
        for (const auto& choice: matrix.getRowGroupIndices(state)) {
            auto row = matrix.getRow(choice);
            uint64_t curCol = 0;
            for (const auto& entry: row) {
                for (; curCol < entry.getColumn(); curCol++) {
                    std::cout << "0 ";
                }
                std::cout << entry.getValue() << " ";
                curCol++;
            }
            for (; curCol < matrix.getColumnCount(); curCol++) {
                std::cout << "0 ";
            }
            std::cout << std::endl << std::endl;
        }
    }
}

template class DARewardProductBuilder<double, storm::models::sparse::StandardRewardModel<double>>;

template class DARewardProductBuilder<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}


#include "DARewardProductBuilder.h"

#include "../modelchecker/CheckTask.h"
#include "../modelchecker/helper/ltl/SparseLTLHelper.h"

#include <gmpxx.h>  //???

#include <numeric>

namespace storm {
namespace transformer {

template<typename ValueType, typename RewardModelType>
std::vector<std::list<uint64_t>> DARewardProductBuilder<ValueType, RewardModelType>::processProductMatrix(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, std::vector<uint64_t> const& stateToMec, storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs, std::vector<uint64_t>& choiceToModelChoice) {

    std::vector<std::list<uint64_t>> mecsToLeavingActions(mecs.size(), std::list<uint64_t>());
    bool isMatrixBuilderEmpty = true;

    // Iterates over all state-action pairs and modifies them for new model
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        uint64_t const mecIndex = stateToMec[state];
        transitionMatrixBuilder.newRowGroup(isMatrixBuilderEmpty ? 0 : transitionMatrixBuilder.getLastRow() + 1);

        for (auto const& choice: transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);
            uint64_t nextRow = isMatrixBuilderEmpty ? 0 : transitionMatrixBuilder.getLastRow() + 1;

            if (mecIndex == InvalidIndex) {
                modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size(), isMatrixBuilderEmpty);
                choiceToModelChoice[transitionMatrixBuilder.getLastRow()] = choice;
                isMatrixBuilderEmpty = false;
                continue;
            }

            if (mecs[mecIndex].containsChoice(state, choice)) {
                // The choice is part of an MEC and therefore cannot leave it
                for (auto& entry: row) {
                    transitionMatrixBuilder.addNextValue(nextRow, entry.getColumn(), entry.getValue());
                }
                choiceToModelChoice[transitionMatrixBuilder.getLastRow()] = choice;
                isMatrixBuilderEmpty = false;
            }   else {
                mecsToLeavingActions[mecIndex].push_back(choice);
            }
        }
    }

    return mecsToLeavingActions;
}


template<typename ValueType, typename RewardModelType>
std::list<uint64_t> DARewardProductBuilder<ValueType, RewardModelType>::addRepresentativeStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, std::vector<uint64_t> const& stateToMec, storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs, std::list<storage::MaximalEndComponent>& accEcs, std::vector<std::list<uint64_t>> const& mecsToLeavingActions) {

    std::list<uint64_t> reachingAccEcChoices;
    storm::storage::BitVector isMECaccepting(mecs.size(), false);
    // add MEC matrix to main matrix
    std::vector<std::list<uint64_t> > mecToMacStates(mecs.size(), std::list<uint64_t>());
    // represents the first state of each MAC so we can later add an action that leads there a.s.
    uint64_t representativeColumn = transitionMatrix.getColumnCount() + mecs.size();

    // determining to which MEC the end component belongs to
    for (auto it=accEcs.begin(); it!=accEcs.end(); ++it) {
        uint64_t representativeEcState = it->begin()->first;
        auto const mec_index = stateToMec[representativeEcState];
        STORM_LOG_ASSERT(mec_index < mecs.size(), "MEC index out of range.");
        std::cout << std::distance(accEcs.begin(), it) << ", " << mec_index << ", " << mecs[mec_index].size() << ", " << it->size() << std::endl;
        // check if the entire MEC is accepting
        if (mecs[mec_index].size() == it->size()) {
            mecToMacStates[mec_index].push_back(representativeEcState);
            isMECaccepting.set(mec_index, true);
            // we do not need to add the accepting ec since the entire MEC is accepting
            it = accEcs.erase(it);
            --it;
        } else {
            mecToMacStates[mec_index].push_back(representativeColumn);
            representativeColumn += it->size();
        }
    }

    // add states representing MECs and the enabled actions that leave them
    for (uint64_t mecIndex = 0; mecIndex < mecs.size(); ++mecIndex) {
        transitionMatrixBuilder.newRowGroup(transitionMatrixBuilder.getLastRow() + 1);

        for (auto const& choice: mecsToLeavingActions[mecIndex]) {
            auto row = transitionMatrix.getRow(choice);
            modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size());
        }

        // actions that reach copies of accepting end components
        for (auto const& macState: mecToMacStates[mecIndex]) {
            uint64_t choice = transitionMatrixBuilder.getLastRow() + 1;
            reachingAccEcChoices.push_back(choice);
            transitionMatrixBuilder.addNextValue(choice, macState, 1);
        }

        // action to representative state of MEC
        if (!isMECaccepting[mecIndex]) {
            auto representativeState = mecs[mecIndex].begin()->first;
            transitionMatrixBuilder.addNextValue(transitionMatrixBuilder.getLastRow() + 1, representativeState, 1);
        }
    }

    if (!accEcs.size()) {
        transitionMatrixBuilder.addNextValue(transitionMatrixBuilder.getLastRow() + 1, transitionMatrix.getColumnCount() + mecs.size() - 1, 1);
    }

    return reachingAccEcChoices;
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::addMACStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder, uint64_t numberMecs, std::list<storage::MaximalEndComponent> accEcs) {
    uint64_t column_offset = transitionMatrix.getColumnCount() + numberMecs;
    // add MACs to matrix
    for (auto const& mac: accEcs) {
        auto stateSet = mac.getStateSet();

        //states and choices from original matrix
        for (auto const& [state, choices]: mac) {
            transitionMatrixBuilder.newRowGroup(transitionMatrixBuilder.getLastRow() + 1);

            for (auto const& choice: choices) {
                auto row = transitionMatrix.getRow(choice);
                uint64_t choiceModified = transitionMatrixBuilder.getLastRow() + 1;

                for (auto& entry: row) {
                    uint64_t nextState = std::distance(stateSet.begin(), stateSet.find(entry.getColumn()));
                    transitionMatrixBuilder.addNextValue(choiceModified, nextState + column_offset, entry.getValue());
                }
            }
        }
        column_offset += stateSet.size();
    }
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::computeConversionsFromModel(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::MaximalEndComponentDecomposition<ValueType> const& mecs, std::list<storage::MaximalEndComponent> const& accEcs, std::vector<uint64_t>& stateToModelState, std::vector<uint64_t>& choiceToModelChoice) {
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        stateToModelState[state] = product.getModelState(state);
    }

    uint64_t stateCounter = transitionMatrix.getRowGroupCount() + mecs.size();
    uint64_t choiceCounter = transitionMatrix.getRowCount() + mecs.size() + accEcs.size();
    // add MACs to matrix
    for (auto const& mac: accEcs) {
        //states and choices from original matrix
        for (auto const& [state, choices]: mac) {
            auto modelState = product.getModelState(state);
            stateToModelState[stateCounter++] = modelState;
            for (auto const& choice: choices) {
                uint64_t rowOffset = choice - transitionMatrix.getRowGroupIndices()[state];
                auto modelChoice = originalModel.getTransitionMatrix().getRowGroupIndices()[modelState] + rowOffset;
                choiceToModelChoice[choiceCounter++] = modelChoice;
            }
        }
    }
}

template<typename ValueType, typename RewardModelType>
std::shared_ptr<DARewardProduct<ValueType>> DARewardProductBuilder<ValueType, RewardModelType>::build() {
    auto transitionMatrix = productModel.getTransitionMatrix();
    auto backwardTransitions = productModel.getBackwardTransitions();

    // MatrixBuilder to build the transition matrix for the demerged MDP
    storage::SparseMatrixBuilder<ValueType> transitionMatrixBuilder(0, 0, 0, false, true, 0);
    //Compute MECs for entire MDP
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(transitionMatrix, backwardTransitions);
    STORM_LOG_INFO(mecs.statistics(transitionMatrix.getRowGroupCount()));
    // Compute accepting end components for the product model

    std::list<storage::MaximalEndComponent> accEcs = computeAcceptingECs(acceptance, transitionMatrix, backwardTransitions);

    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), InvalidIndex);
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        for (auto const& [state, _] : mec) {
            stateToMec[state] = mec_counter;
        }
        ++mec_counter;
    }

    uint64_t numberOfChoicesAccEcs = 0, numberOfStatesAccEcs = 0;
    for (auto const& ec: accEcs) {
        for (auto const& [_, choices] : ec) {
            numberOfChoicesAccEcs += choices.size();
        }
        numberOfStatesAccEcs += ec.size();
    }

    STORM_LOG_INFO("Found " << accEcs.size() << " accepting end components with a total number of " << numberOfStatesAccEcs << " states and " << numberOfChoicesAccEcs << " choices." );
    std::vector<uint64_t> stateToModelState(numberOfStatesAccEcs + transitionMatrix.getRowGroupCount() + mecs.size(), InvalidIndex);
    std::vector<uint64_t> choiceToModelChoice(numberOfChoicesAccEcs + transitionMatrix.getRowCount() + mecs.size() + accEcs.size(), InvalidIndex);

    // Processes the stated from the product model for the modified one
    auto mecsToLeavingActions = processProductMatrix(transitionMatrix, transitionMatrixBuilder, stateToMec, mecs, choiceToModelChoice);
    // Adds states representing MECs to the matrix builder
    auto reachingAccEcChoices = addRepresentativeStates(transitionMatrix, transitionMatrixBuilder, stateToMec, mecs, accEcs, mecsToLeavingActions);

    addMACStates(transitionMatrix, transitionMatrixBuilder, mecs.size(), accEcs);

    auto modifiedTransitionMatrix = transitionMatrixBuilder.build();
    STORM_LOG_ASSERT(modifiedTransitionMatrix.isProbabilistic(), "The resulting transition matrix is not probabilistic");
    STORM_LOG_INFO("The modified model has " << modifiedTransitionMatrix.getRowGroupCount() << " states and " << modifiedTransitionMatrix.getRowCount() << " choices.");
    auto initialStates = liftInitialStates(stateToMec, modifiedTransitionMatrix.getRowGroupCount());
    computeConversionsFromModel(transitionMatrix, mecs, accEcs, stateToModelState, choiceToModelChoice);

    return std::make_shared<DARewardProduct<ValueType>>(modifiedTransitionMatrix, stateToModelState, choiceToModelChoice, reachingAccEcChoices, initialStates);
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, std::vector<uint64_t> const& stateToMec, Row row, uint64_t numMecs, bool builderEmpty) {
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

    removeSubsetEndComponents(macs);
    return macs;
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::removeSubsetEndComponents(std::list<storage::MaximalEndComponent>& endComponents) {
    for (auto it = endComponents.begin(); it != endComponents.end(); ++it) {
        const auto& stateSetA = it->getStateSet();

        auto nextIt = it;
        ++nextIt;

        while (nextIt != endComponents.end()) {
            const auto& stateSetB = nextIt->getStateSet();

            if (std::includes(stateSetA.begin(), stateSetA.end(), stateSetB.begin(), stateSetB.end())) {
                // stateSetB is a subset of stateSetA
                nextIt = endComponents.erase(nextIt);
            } else if (std::includes(stateSetB.begin(), stateSetB.end(), stateSetA.begin(), stateSetA.end())) {
                // stateSetA is a subset of stateSetB
                it = endComponents.erase(it);
                --it;
                break;
            } else {
                ++nextIt;
            }
        }
    }
}


template<typename ValueType, typename RewardModelType>
storm::storage::BitVector DARewardProductBuilder<ValueType, RewardModelType>::liftInitialStates(std::vector<uint64_t> const& stateToMec, uint64_t const numberOfStates) {
    storm::storage::BitVector initialStatesModified(numberOfStates, false);

    for (uint64_t prodInitState: initialStatesProduct) {
        uint64_t mecIndex = stateToMec[prodInitState];
        if (mecIndex == InvalidIndex) {
            initialStatesModified.set(prodInitState, true);
        } else {
            initialStatesModified.set(stateToMec.size() + mecIndex, true);
        }
    }
    return initialStatesModified;
}

template class DARewardProductBuilder<double, storm::models::sparse::StandardRewardModel<double>>;

template class DARewardProductBuilder<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}
}


#include "DARewardProductBuilder.h"

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/automata/AcceptanceConditionSynthesizer.h"
#include "storm/modelchecker/CheckTask.h"
#include "storm/modelchecker/helper/ltl/SparseLTLHelper.h"

#include <numeric>

namespace storm {
namespace transformer {

template<typename ValueType>
void printCompactTransitionMatrix(const storm::storage::SparseMatrix<ValueType>& transitionMatrix) {
    std::cout << "Transition Matrix:" << std::endl;
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); ++state) {
        std::cout << "  State " << state << ":" << std::endl;

        for (auto const& choice : transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);

            std::cout << "    Choice " << choice << ": ";
            uint64_t currentColumn = 0;

            for (auto const& entry : row) {
                // Fülle Lücken mit Nullen bis zur nächsten belegten Spalte
                for (; currentColumn < entry.getColumn(); ++currentColumn) {
                    std::cout << "0 ";
                }
                // Gebe Wert aus
                std::cout << entry.getValue() << " ";
                ++currentColumn;
            }

            for (; currentColumn < transitionMatrix.getColumnCount(); ++currentColumn) {
                std::cout << "0 ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

template<typename ValueType, typename RewardModelType>
std::vector<std::list<uint64_t>> DARewardProductBuilder<ValueType, RewardModelType>::processProductMatrix(
    storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
    MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions) {
    std::vector<std::list<uint64_t>> mecsToLeavingActions(mecDecompositionInfo.mecs.size(), std::list<uint64_t>());
    bool isMatrixBuilderEmpty = true;

    for (uint64_t state = 0, newState = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        if (mecDecompositionInfo.stateToMec[state] != InvalidIndex)
            continue;
        conversions.modelStateToState[state] = newState;
        conversions.stateToModelState[newState] = productToModelState[state];
        newState++;
    }

    // Iterates over all state-action pairs and modifies them for new model
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        uint64_t const mecIndex = mecDecompositionInfo.stateToMec[state];
        if (mecIndex == InvalidIndex) {
            transitionMatrixBuilder.newRowGroup(isMatrixBuilderEmpty ? 0 : transitionMatrixBuilder.getLastRow() + 1);
        }

        for (auto const& choice : transitionMatrix.getRowGroupIndices(state)) {
            if (mecIndex == InvalidIndex) {
                modifyStateActionPair(transitionMatrixBuilder, transitionMatrix.getRow(choice), mecDecompositionInfo, conversions, isMatrixBuilderEmpty);
                conversions.choiceToModelChoice[transitionMatrixBuilder.getLastRow()] = choice;
                isMatrixBuilderEmpty = false;
                continue;
            }

            if (!mecDecompositionInfo.mecs[mecIndex].containsChoice(state, choice)) {
                mecsToLeavingActions[mecIndex].push_back(choice);
            }
        }
    }

    return mecsToLeavingActions;
}

template<typename ValueType, typename RewardModelType>
std::vector<std::list<uint64_t>> DARewardProductBuilder<ValueType, RewardModelType>::addRepresentativeStates(
    storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
    std::vector<std::list<storage::MaximalEndComponent>> const& accEcs, std::vector<std::list<uint64_t>> const& mecsToLeavingActions,
    MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions) {
    std::vector reachingAccEcChoices(accEcs.size(), std::list<uint64_t>());
    uint64_t numMECs = mecsToLeavingActions.size();
    // add MEC matrix to main matrix
    std::vector mecToMacStates(numMECs, std::vector(accEcs.size(), std::list<uint64_t>()));
    // represents the first state of each MAC so we can later add an action that leads there a.s.
    uint64_t representativeColumn = transitionMatrix.getRowGroupCount() - mecDecompositionInfo.numStatesInMecs + numMECs;

    for (int i = 0; i < accEcs.size(); i++) {
        // determining to which MEC the end component belongs to
        for (auto const& ec : accEcs[i]) {
            uint64_t representativeEcState = ec.begin()->first;
            auto const mec_index = mecDecompositionInfo.stateToMec[representativeEcState];
            STORM_LOG_ASSERT(mec_index < numMECs, "MEC index out of range.");
            mecToMacStates[mec_index][i].push_back(representativeColumn);
            representativeColumn += ec.size();
        }
    }

    bool isMatrixBuilderEmpty = mecDecompositionInfo.numStatesInMecs == mecDecompositionInfo.stateToMec.size();
    // add states representing MECs and the enabled actions that leave them
    for (uint64_t mecIndex = 0; mecIndex < numMECs; ++mecIndex) {
        transitionMatrixBuilder.newRowGroup(isMatrixBuilderEmpty ? 0 : transitionMatrixBuilder.getLastRow() + 1);

        // actions that might reach states not in some end component
        for (auto const& choice : mecsToLeavingActions[mecIndex]) {
            modifyStateActionPair(transitionMatrixBuilder, transitionMatrix.getRow(choice), mecDecompositionInfo, conversions);
        }

        for (int i = 0; i < accEcs.size(); i++) {
            // actions that reach copies of accepting end components
            for (auto const& macState : mecToMacStates[mecIndex][i]) {
                uint64_t choice = isMatrixBuilderEmpty ? 0 : transitionMatrixBuilder.getLastRow() + 1;
                isMatrixBuilderEmpty = false;
                reachingAccEcChoices[i].push_back(choice);
                transitionMatrixBuilder.addNextValue(choice, macState, 1);
            }
        }
    }

    return reachingAccEcChoices;
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::addMACStates(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                      storage::SparseMatrixBuilder<ValueType>& transitionMatrixBuilder,
                                                                      std::list<storage::MaximalEndComponent> const& accEcs,
                                                                      MecDecompositionInfo mecDecompositionInfo) {
    uint64_t column_offset = transitionMatrixBuilder.getCurrentRowGroupCount();
    // add MACs to matrix
    for (auto const& mac : accEcs) {
        auto stateSet = mac.getStateSet();

        // states and choices from original matrix
        for (auto const& state : stateSet) {
            transitionMatrixBuilder.newRowGroup(transitionMatrixBuilder.getLastRow() + 1);

            for (auto const& choice : mac.getChoicesForState(state)) {
                auto row = transitionMatrix.getRow(choice);
                uint64_t choiceModified = transitionMatrixBuilder.getLastRow() + 1;
                for (auto& entry : row) {
                    uint64_t nextState = std::distance(stateSet.begin(), stateSet.find(entry.getColumn()));
                    transitionMatrixBuilder.addNextValue(choiceModified, nextState + column_offset, entry.getValue());
                }
            }
        }
        column_offset += stateSet.size();
    }
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::computeConversionsFromModel(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                                     std::vector<std::list<storage::MaximalEndComponent>> const& accEcs,
                                                                                     MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions,
                                                                                     uint64_t stateCounter, uint64_t choiceCounter) {
    std::list<storage::MaximalEndComponent> flattenedEcs;
    for (auto ecs : accEcs) {
        flattenedEcs.splice(flattenedEcs.end(), ecs);
    }
    std::cout << flattenedEcs.size() << std::endl;
    // add MACs to matrix
    for (auto const& mac : flattenedEcs) {
        // states and choices from original matrix
        for (auto const& state : mac.getStateSet()) {
            auto modelState = productToModelState[state];
            conversions.stateToModelState[stateCounter++] = modelState;
            for (auto const& choice : mac.getChoicesForState(state)) {
                uint64_t rowOffset = choice - transitionMatrix.getRowGroupIndices()[state];
                auto modelChoice = originalModel.getTransitionMatrix().getRowGroupIndices()[modelState] + rowOffset;
                conversions.choiceToModelChoice[choiceCounter++] = modelChoice;
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
    // Compute MECs for entire MDP
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(transitionMatrix, backwardTransitions);
    STORM_LOG_INFO(mecs.statistics(transitionMatrix.getRowGroupCount()));
    // Compute accepting end components for the product model
    uint64_t numberAcceptanceConditionCombinations = pow(2, acceptanceConditions.size());

    auto objectivesCombinations = storm::automata::AcceptanceConditionSynthesizer::getAllCombinations(acceptanceConditions);
    std::vector<std::list<storage::MaximalEndComponent>> accEcs(numberAcceptanceConditionCombinations);

    // computes the accepting end components for all combinations of the LTL formula
    for (uint64_t i = 0; i < numberAcceptanceConditionCombinations; ++i) {
        accEcs[i] = computeAcceptingECs(*objectivesCombinations[i], transitionMatrix, backwardTransitions);
        STORM_LOG_INFO("Found " << accEcs[i].size() << " accepting end components for combination " << i);
    }

    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), InvalidIndex);
    uint64_t numStatesInMec = 0, numChoicesInMec = 0;
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        for (auto const& [state, choices] : mec) {
            stateToMec[state] = mec_counter;
            numChoicesInMec += choices.size();
        }
        numStatesInMec += mec.size();
        ++mec_counter;
    }
    MecDecompositionInfo mecDecompositionInfo(mecs, numStatesInMec, numChoicesInMec, stateToMec);
    uint64_t numberOfChoicesAccEcs = 0, numberOfStatesAccEcs = 0, numberOfAccEcs = 0;
    for (uint64_t i = 0; i < numberAcceptanceConditionCombinations; i++) {
        for (auto const& ec : accEcs[i]) {
            for (auto const& [_, choices] : ec) {
                numberOfChoicesAccEcs += choices.size();
            }
            numberOfStatesAccEcs += ec.size();
            numberOfAccEcs++;
        }
    }

    STORM_LOG_INFO("Found " << accEcs.size() << " accepting end components with a total number of " << numberOfStatesAccEcs << " states and "
                            << numberOfChoicesAccEcs << " choices.");
    const uint64_t totalNumberOfStates = numberOfStatesAccEcs + transitionMatrix.getRowGroupCount() - numStatesInMec + mecs.size();
    const uint64_t totalNumberOfChoices = numberOfChoicesAccEcs + transitionMatrix.getRowCount() - numChoicesInMec + numberOfAccEcs;
    Conversions conversions;
    conversions.stateToModelState.resize(totalNumberOfStates, InvalidIndex);
    conversions.choiceToModelChoice.resize(totalNumberOfChoices, InvalidIndex);
    conversions.modelStateToState.resize(transitionMatrix.getRowGroupCount(), InvalidIndex);

    // Processes the stated from the product model for the modified one
    auto mecsToLeavingActions = processProductMatrix(transitionMatrix, transitionMatrixBuilder, mecDecompositionInfo, conversions);
    // printCompactTransitionMatrix(transitionMatrixBuilder.build());
    auto reachingAccEcChoices =
        addRepresentativeStates(transitionMatrix, transitionMatrixBuilder, accEcs, mecsToLeavingActions, mecDecompositionInfo, conversions);

    for (uint64_t i = 0; i < numberAcceptanceConditionCombinations; i++) {
        addMACStates(transitionMatrix, transitionMatrixBuilder, accEcs[i], mecDecompositionInfo);
        STORM_LOG_INFO("Added " << accEcs[i].size() << " MACs to the modified model.");
    }

    auto modifiedTransitionMatrix = transitionMatrixBuilder.build();
    STORM_LOG_ASSERT(modifiedTransitionMatrix.isProbabilistic(), "The resulting transition matrix is not probabilistic");
    STORM_LOG_ASSERT(modifiedTransitionMatrix.getRowGroupCount() == totalNumberOfStates, "The resulting transition matrix has the wrong number of states");
    STORM_LOG_ASSERT(modifiedTransitionMatrix.getRowCount() == totalNumberOfChoices, "The resulting transition matrix has the wrong number of choices");
    STORM_LOG_INFO("The modified model has " << modifiedTransitionMatrix.getRowGroupCount() << " states and " << modifiedTransitionMatrix.getRowCount()
                                             << " choices.");
    auto initialStates = liftInitialStates(mecDecompositionInfo, totalNumberOfStates, productModel.getInitialStates());

    uint64_t stateCounter = totalNumberOfStates - numberOfStatesAccEcs;
    uint64_t choiceCounter = totalNumberOfChoices - numberOfChoicesAccEcs;
    computeConversionsFromModel(transitionMatrix, accEcs, mecDecompositionInfo, conversions, stateCounter, choiceCounter);

    return std::make_shared<DARewardProduct<ValueType>>(modifiedTransitionMatrix, conversions.stateToModelState, conversions.choiceToModelChoice,
                                                        reachingAccEcChoices, initialStates);
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, Row row,
                                                                               MecDecompositionInfo const& mecDecompositionInfo, Conversions& conversions,
                                                                               bool builderEmpty) {
    std::vector<ValueType> mecValues(mecDecompositionInfo.mecs.size(), 0);
    uint64_t choice = builderEmpty ? 0 : builder.getLastRow() + 1;
    const uint64_t numStates = mecDecompositionInfo.stateToMec.size() - mecDecompositionInfo.numStatesInMecs;

    for (auto& entry : row) {
        uint64_t mecIndexNextState = mecDecompositionInfo.stateToMec[entry.getColumn()];
        if (mecIndexNextState == std::numeric_limits<uint64_t>::max()) {
            uint64_t nextState = conversions.modelStateToState[entry.getColumn()];
            builder.addNextValue(choice, nextState, entry.getValue());
        } else {
            mecValues[mecIndexNextState] += entry.getValue();
        }
    }

    for (uint64_t i = 0; auto& entry : mecValues) {
        if (entry) {
            builder.addNextValue(choice, numStates + i, entry);
        }
        i++;
    }
}

template<typename ValueType, typename RewardModelType>
std::list<storage::MaximalEndComponent> DARewardProductBuilder<ValueType, RewardModelType>::computeAcceptingECs(
    automata::AcceptanceCondition const& acceptance, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions) {
    // list of maximum accepting end components in the product MDP
    std::list<storage::MaximalEndComponent> macs;
    std::vector<std::vector<automata::AcceptanceCondition::acceptance_expr::ptr>> dnf = acceptance.extractFromDNF();

    for (auto const& conjunction : dnf) {
        // get the states of the mdp that are on a MEC and don't violate Fins of the conjunction
        storm::storage::BitVector allowed(transitionMatrix.getRowGroupCount(), true);  // maybe late start with mecStates

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
storm::storage::BitVector DARewardProductBuilder<ValueType, RewardModelType>::liftInitialStates(MecDecompositionInfo const& mecDecompositionInfo,
                                                                                                uint64_t numberOfStates,
                                                                                                storm::storage::BitVector const& initialStates) const {
    storm::storage::BitVector initialStatesModified(numberOfStates);
    uint64_t offset = mecDecompositionInfo.stateToMec.size() - mecDecompositionInfo.numStatesInMecs;

    for (uint64_t prodInitState : initialStates) {
        uint64_t mecIndex = mecDecompositionInfo.stateToMec[prodInitState];
        if (mecIndex == InvalidIndex) {
            initialStatesModified.set(prodInitState, true);
        } else {
            initialStatesModified.set(offset + mecIndex, true);
        }
    }
    return initialStatesModified;
}

template class DARewardProductBuilder<double, storm::models::sparse::StandardRewardModel<double>>;

template class DARewardProductBuilder<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;

}  // namespace transformer
}  // namespace storm

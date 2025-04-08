#include "DARewardProductBuilder.h"

#include "../modelchecker/CheckTask.h"
#include "../modelchecker/helper/ltl/SparseLTLHelper.h"
#include "../models/sparse/MarkovAutomaton.h"
#include "../storage/prism/RewardModel.h"
#include "storm/automata/DeterministicAutomaton.h"
#include "storm/automata/LTL2DeterministicAutomaton.h"
#include "storm/logic/ExtractMaximalStateFormulasVisitor.h"

#include "storm/exceptions/InvalidPropertyException.h"
#include <gmpxx.h>  //???

#define DEBUG_VAR(x) std::cout << #x << " = " << x << std::endl;

namespace storm {
namespace transformer {

template <typename T>
void printMatrix(storm::storage::SparseMatrix<T> matrix);

template<typename ValueType, typename RewardModelType>
std::shared_ptr<DAProduct<models::sparse::Mdp<ValueType, RewardModelType>>>  DARewardProductBuilder<ValueType, RewardModelType>::buildProductMDP() {
    storm::logic::ExtractMaximalStateFormulasVisitor::ApToFormulaMap extracted;
    std::shared_ptr<storm::logic::Formula> ltlFormula = storm::logic::ExtractMaximalStateFormulasVisitor::extract(formula, extracted);
    STORM_LOG_ASSERT(ltlFormula->isPathFormula(), "Unexpected formula type.");

    // Compute Satisfaction sets for the APs (which represent the state-subformulae
    auto apSatSets = storm::modelchecker::helper::SparseLTLHelper<ValueType, true>::computeApSets(extracted, formulaChecker);

    STORM_LOG_INFO("Resulting LTL path formula: " << ltlFormula->toString());

    std::shared_ptr<storm::automata::DeterministicAutomaton> da = storm::automata::LTL2DeterministicAutomaton::ltl2daSpot(*ltlFormula, true);

    return buildDAProduct(*da, apSatSets);
}


template<typename ValueType, typename RewardModelType>
std::shared_ptr<DAProduct<models::sparse::Mdp<ValueType, RewardModelType>>> DARewardProductBuilder<ValueType, RewardModelType>::buildDAProduct(storm::automata::DeterministicAutomaton const& da, std::map<std::string, storm::storage::BitVector>& apSatSets) {
    const storm::automata::APSet& apSet = da.getAPSet();

    std::vector<storm::storage::BitVector> statesForAP;
    for (const std::string& ap : apSet.getAPs()) {
        auto it = apSatSets.find(ap);
        STORM_LOG_THROW(it != apSatSets.end(), storm::exceptions::InvalidOperationException,
                        "Deterministic automaton has AP " << ap << ", does not appear in formula");

        statesForAP.push_back(std::move(it->second));
    }

    storm::storage::BitVector statesOfInterest(model.getNumberOfStates(), true);
    transformer::DAProductBuilder productBuilder(da, statesForAP);

    auto product = productBuilder.build<Mdp>(model.getTransitionMatrix(), statesOfInterest);

    STORM_LOG_INFO("Product MDP-DA has "
                   << product->getProductModel().getNumberOfStates() << " states and " << product->getProductModel().getNumberOfTransitions()
                   << " transitions.");

    return product;
}


template<typename ValueType, typename RewardModelType>
std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, storage::SparseMatrix<ValueType>> DARewardProductBuilder<ValueType, RewardModelType>::buildTransitionMatrix(automata::AcceptanceCondition const& acceptance, storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::storage::SparseMatrix<ValueType> const& backwardTransitions) {
    STORM_LOG_INFO("Building transition matrix...");

    // map every state/action in the modified model to the original state/action in the product model
    std::vector<uint64_t> stateToProductState;
    std::vector<uint64_t> actionToOriginalAction;
    // MatrixBuilder to build the transition matrix for the demerged MDP
    storage::SparseMatrixBuilder<ValueType> transitionMatrixBuilder(0, 0, 0, false, true, 0);
    auto lastRow = [&transitionMatrixBuilder](){ return transitionMatrixBuilder.getLastRow();};

    //Compute MECs for entire MDP
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(transitionMatrix, backwardTransitions);
    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    uint64_t InvalidMECIndex = std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), InvalidMECIndex);
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        std::cout << "MEC contains:" << std::endl;
        for (auto const& [state, _] : mec) {
            stateToMec[state] = mec_counter;
        }
        ++mec_counter;
    }

    // compute the accepting end components
    std::list<storage::MaximalEndComponent> accEcs = computeAcceptingECs(acceptance, transitionMatrix, backwardTransitions);
    uint64_t totalNumberStates = transitionMatrix.getRowGroupCount() + mecs.size();
    for (auto& mac: accEcs) {
        totalNumberStates += mac.size();
    }
    stateToProductState.reserve(totalNumberStates);

    // A mapping of each mec to the actions that are enabled in a state in it but leave the mec with some probability
    std::vector<std::list<uint64_t>> mecsToLeavingActions(mecs.size(), std::list<uint64_t>());
    bool isMatrixBuilderEmpty = true;

    // Iterates over all state-action pairs and modifies them for new model
    for (uint64_t state = 0; state < transitionMatrix.getRowGroupCount(); state++) {
        uint64_t const mecIndex = stateToMec[state];
        uint64_t startRowGroup = isMatrixBuilderEmpty ? 0 : lastRow() + 1;
        transitionMatrixBuilder.newRowGroup(startRowGroup);

        for (auto const& choice: transitionMatrix.getRowGroupIndices(state)) {
            auto row = transitionMatrix.getRow(choice);
            uint64_t nextRow = isMatrixBuilderEmpty ? 0 : lastRow() + 1;

            DEBUG_VAR(choice);

            if (mecIndex == InvalidMECIndex) {
                modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size(), isMatrixBuilderEmpty);
            }

            bool choiceIsInMec = (mecIndex != InvalidMECIndex) && mecs[mecIndex].containsChoice(state, choice);

            // The choice is part of an MEC and therefore cannot leave it
            if (choiceIsInMec) {
                //stateToMec[state] = nextRow;
                for (auto& entry: row) {
                    transitionMatrixBuilder.addNextValue(nextRow, entry.getColumn(), entry.getValue());
                }
            }

            if (mecIndex != InvalidMECIndex && !choiceIsInMec) {
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

    // determining to which MEC the end component belong to
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

        for (auto const& choice: mecsToLeavingActions[mecIndex]) {
            auto row = transitionMatrix.getRow(choice);
            modifyStateActionPair(transitionMatrixBuilder, stateToMec, row, mecs.size());
        }

        for (auto const& macState: mecToMacStates[mecIndex]) {
            transitionMatrixBuilder.addNextValue(lastRow() + 1, macState, 1);
        }

        // action to representative state of MEC
        auto representativeState = mecs[mecIndex].begin()->first;
        transitionMatrixBuilder.addNextValue(lastRow() + 1, representativeState, 1);
    }

    uint64_t column_offset = transitionMatrix.getColumnCount() + mecs.size();
    // add MACs to matrix
    for (auto const& mac: accEcs) {
        auto stateSet = mac.getStateSet();

        //states and choices from original matrix
        for (auto const& [state, choices]: mac) {
            transitionMatrixBuilder.newRowGroup(lastRow() + 1);

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

    return std::make_tuple(stateToProductState, actionToOriginalAction, transitionMatrixBuilder.build());
}

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::modifyStateActionPair(storage::SparseMatrixBuilder<ValueType>& builder, std::vector<uint64_t>& stateToMec, typename storage::SparseMatrix<ValueType>::const_rows row, uint64_t numMecs, bool builderEmpty) {
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

    //ensure matrix is probabilistic
    //STORM_LOG_ASSERT(sumEntries == 1, "Sum of entries in row do not add up to one");
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

template<typename ValueType, typename RewardModelType>
void DARewardProductBuilder<ValueType, RewardModelType>::build() {
    std::cout << "Building the demerged matrix..." << std::endl;

    auto product = buildProductMDP();
    std::cout << "The product matrix..." << std::endl;
    printMatrix(product->getProductModel().getTransitionMatrix());

    auto result = buildTransitionMatrix(*product->getAcceptance(), product->getProductModel().getTransitionMatrix(), product->getProductModel().getBackwardTransitions());

    auto transitionMatrix = std::get<2>(result);
    std::cout << "Resulting transition matrix:" << std::endl;
    printMatrix(transitionMatrix);
    STORM_LOG_ASSERT(transitionMatrix.isProbabilistic(), "The resulting transition matrix is not probabilistic");
    std::cout << "Done building the demerged matrix." << std::endl;
}


template class DARewardProductBuilder<double, storm::models::sparse::StandardRewardModel<double>>;

template class DARewardProductBuilder<RationalNumber, storm::models::sparse::StandardRewardModel<RationalNumber>>;



#ifdef DA_READY


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




#endif
}
}


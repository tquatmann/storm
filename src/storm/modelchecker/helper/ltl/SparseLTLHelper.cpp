#include "SparseLTLHelper.h"

#include "storm/automata/AcceptanceConditionSynthesizer.h"
#include "storm/automata/Automaton.h"
#include "storm/automata/LTL2Automaton.h"

#include "storm/environment/modelchecker/ModelCheckerEnvironment.h"

#include "storm/logic/ExtractMaximalStateFormulasVisitor.h"

#include "storm/modelchecker/prctl/helper/SparseDtmcPrctlHelper.h"
#include "storm/modelchecker/prctl/helper/SparseMdpPrctlHelper.h"
#include "storm/modelchecker/propositional/SparsePropositionalModelChecker.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"

#include "storm/solver/SolveGoal.h"
#include "storm/storage/MaximalEndComponentDecomposition.h"
#include "storm/storage/SchedulerChoice.h"
#include "storm/storage/StronglyConnectedComponentDecomposition.h"

#include "storm/exceptions/InvalidPropertyException.h"
#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"

namespace storm {
namespace modelchecker {
namespace helper {

namespace detail {

/*!
 * Builds an automaton for the given formula.
 * @tparam Automaton
 * @param formula
 * @param env
 * @return
 */
template<typename Automaton>
typename Automaton::ptr buildAutomatonFromFormula(Environment const& env, storm::logic::Formula const& ltlFormula, bool forceDnfAcceptance) {
    typename Automaton::ptr automaton;
    if (env.modelchecker().isLtl2AutToolSet()) {
        automaton = storm::automata::LTL2Automaton::ltl2AutExternalTool<Automaton>(ltlFormula, env.modelchecker().getLtl2AutTool());
        if (auto accExprPtr = automaton->getAcceptance()->getAcceptanceExpression();
            forceDnfAcceptance && !storm::automata::isDisjunctiveNormalForm(accExprPtr)) {
            STORM_LOG_WARN("Converting acceptance condition "
                           << *accExprPtr
                           << " into disjunctive normal form. This might blow up the size of the acceptance condition. Check if the LTL2Aut "
                              "tool can be configured to produce DNFs directly.");
            automaton->getAcceptance()->setAcceptanceExpression(storm::automata::toDisjunctiveNormalForm(accExprPtr));
        }
    } else {
        if constexpr (std::is_same_v<Automaton, storm::automata::DeterministicAutomaton>) {
            automaton = storm::automata::LTL2Automaton::ltl2AutSpot(ltlFormula, forceDnfAcceptance);
            STORM_LOG_THROW(!forceDnfAcceptance || storm::automata::isDisjunctiveNormalForm(automaton->getAcceptance()->getAcceptanceExpression()),
                            storm::exceptions::UnexpectedException, "Spot produced an automaton whose acceptance condition is not in disjunctive normal form.");
        } else {
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "The internal engine (Spot) can only generate deterministic automata.");
        }
    }
    STORM_LOG_INFO(Automaton::type << " for LTL formula has " << automaton->getNumberOfStates() << " states, " << automaton->getAPSet().size()
                                   << " atomic propositions and " << *automaton->getAcceptance()->getAcceptanceExpression() << " as acceptance condition.\n");
    return automaton;
}

/*!
 * Returns for each AP the set of states that satisfy it, using the order given by the provided apSet
 * @param indexMapping if given, the resulting state sets are with respect to the indices given in the mapping, i.e., the returned BitVectors will have
 * indexMapping->size() entries where the i'th bit matches the original bit for indexMapping->at(i)
 */
std::vector<storm::storage::BitVector> computeStatesForApSet(storm::automata::APSet const& apSet, std::map<std::string, storm::storage::BitVector>&& apSatSets,
                                                             std::optional<std::vector<uint64_t> const> indexMapping = {}) {
    std::vector<storm::storage::BitVector> result;
    for (auto const& ap : apSet.getAPs()) {
        auto it = apSatSets.find(ap);
        STORM_LOG_ASSERT(it != apSatSets.end(), "AP " << ap << " of automaton does not appear in formula");
        if (indexMapping) {
            auto& newSet = result.emplace_back(indexMapping->size(), false);
            for (uint64_t i = 0; i < newSet.size(); ++i) {
                uint64_t oldIndex = (*indexMapping)[i];
                if (it->second.get(oldIndex)) {
                    result.back().set(i);
                }
            }
        } else {
            result.push_back(std::move(it->second));
        }
    }
    return result;
}

RabinObjective extractRabinObjectiveFromProductAcceptance(storm::automata::AcceptanceCondition&& acceptance, storm::models::sparse::StateLabeling& labeling) {
    RabinObjective result;
    std::vector<std::vector<automata::AcceptanceCondition::acceptance_expr::ptr>> dnf = acceptance.extractFromDNF();
    for (auto const& conjunction : dnf) {
        if (std::any_of(conjunction.begin(), conjunction.end(), [](auto const& literal) { return literal->isFALSE(); })) {
            continue;  // ignore conjunction in expression "(X && Y && false && Z) || ..."
        }
        auto& rabinCond = result.emplace_back();
        for (auto const& literal : conjunction) {
            if (literal->isTRUE()) {
                continue;  // ignore true in expression "true && ..."
            } else if (literal->isAtom()) {
                auto atom = literal->getAtom();
                storm::storage::BitVector stateSet;
                if (atom.isNegated()) {
                    stateSet = ~acceptance.getAcceptanceSet(atom.getAcceptanceSet());
                } else {
                    stateSet = std::move(acceptance.getAcceptanceSet(atom.getAcceptanceSet()));
                }
                std::string apName;
                if (auto foundLabel = labeling.findLabel(stateSet); foundLabel.has_value()) {
                    apName = *foundLabel;
                } else {
                    apName = labeling.addUniqueLabel("rabin_obj_ap", std::move(stateSet));
                }
                if (atom.getType() == cpphoafparser::AtomAcceptance::TEMPORAL_FIN) {
                    rabinCond.fin.push_back(apName);
                } else {
                    STORM_LOG_ASSERT(atom.getType() == cpphoafparser::AtomAcceptance::TEMPORAL_INF, "Unknown acceptance atom type.");
                    rabinCond.inf.push_back(apName);
                }
            } else {
                STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected acceptance literal type.");
            }
        }
    }
    return result;
}
}  // namespace detail

template<typename ValueType, bool Nondeterministic>
SparseLTLHelper<ValueType, Nondeterministic>::SparseLTLHelper(storm::storage::SparseMatrix<ValueType> const& transitionMatrix)
    : _transitionMatrix(transitionMatrix) {
    // Intentionally left empty.
}

template<typename ValueType, bool Nondeterministic>
storm::storage::Scheduler<ValueType> SparseLTLHelper<ValueType, Nondeterministic>::SparseLTLHelper::extractScheduler(
    storm::models::sparse::Model<ValueType> const& model) {
    STORM_LOG_ASSERT(this->isProduceSchedulerSet(), "Trying to get the produced optimal choices although no scheduler was requested.");
    STORM_LOG_ASSERT(this->_schedulerHelper.is_initialized(),
                     "Trying to get the produced optimal choices but none were available. Was there a computation call before?");

    return this->_schedulerHelper.get().extractScheduler(model, this->hasRelevantStates());
}

template<typename ValueType, bool Nondeterministic>
std::vector<ValueType> SparseLTLHelper<ValueType, Nondeterministic>::computeLTLProbabilities(Environment const& env, storm::logic::PathFormula const& formula,
                                                                                             CheckFormulaCallback const& formulaChecker) {
    // Replace state-subformulae by atomic propositions (APs)
    storm::logic::ExtractMaximalStateFormulasVisitor::ApToFormulaMap extracted;
    std::shared_ptr<storm::logic::Formula> ltlFormula = storm::logic::ExtractMaximalStateFormulasVisitor::extract(formula, extracted);
    STORM_LOG_ASSERT(ltlFormula->isPathFormula(), "Unexpected formula type.");

    // Compute Satisfaction sets for the APs (which represent the state-subformulae
    auto apSets = computeApSets(extracted, formulaChecker);

    STORM_LOG_INFO("Computing LTL probabilities for formula with " << apSets.size() << " atomic proposition(s).");

    // Compute the resulting LTL probabilities
    return computeLTLProbabilities(env, ltlFormula->asPathFormula(), std::move(apSets));
}

template<typename ValueType, bool Nondeterministic>
std::map<std::string, storm::storage::BitVector> SparseLTLHelper<ValueType, Nondeterministic>::computeApSets(
    std::map<std::string, std::shared_ptr<storm::logic::Formula const>> const& extracted, CheckFormulaCallback const& formulaChecker) {
    std::map<std::string, storm::storage::BitVector> apSets;
    for (auto& p : extracted) {
        STORM_LOG_DEBUG(" Computing satisfaction set for atomic proposition \"" << p.first << "\" <=> " << *p.second << "...");
        apSets[p.first] = formulaChecker(*p.second);
    }
    return apSets;
}

template<typename ValueType, bool Nondeterministic>
storm::storage::BitVector SparseLTLHelper<ValueType, Nondeterministic>::computeAcceptingECs(automata::AcceptanceCondition const& acceptance,
                                                                                            storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                                            storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                                                            typename transformer::DAProduct<productModelType>::ptr product) {
    STORM_LOG_INFO("Computing accepting states for acceptance condition " << *acceptance.getAcceptanceExpression());
    if (acceptance.getAcceptanceExpression()->isTRUE()) {
        STORM_LOG_INFO(" TRUE -> all states accepting (assumes no deadlock in the model)");
        return storm::storage::BitVector(transitionMatrix.getRowGroupCount(), true);
    } else if (acceptance.getAcceptanceExpression()->isFALSE()) {
        STORM_LOG_INFO(" FALSE -> all states rejecting");
        return storm::storage::BitVector(transitionMatrix.getRowGroupCount(), false);
    }

    std::vector<std::vector<automata::AcceptanceCondition::acceptance_expr::ptr>> dnf = acceptance.extractFromDNF();

    storm::storage::BitVector acceptingStates(transitionMatrix.getRowGroupCount(), false);

    std::size_t accMECs = 0;
    std::size_t allMECs = 0;

    // All accepting states will be on a MEC. For efficiency, we compute the MECs of the MDP first to restrict the possible candidates.
    // Compute MECs for the entire MDP
    storm::storage::MaximalEndComponentDecomposition<ValueType> mecs(transitionMatrix, backwardTransitions);
    // Maps every state to the MEC it is in, or to InvalidMecIndex if it does not belong to any MEC.
    std::vector<uint64_t> stateToMec(transitionMatrix.getRowGroupCount(), std::numeric_limits<uint64_t>::max());
    // Contains states that are on a mec
    storm::storage::BitVector mecStates(transitionMatrix.getRowGroupCount(), false);
    for (uint64_t mec_counter = 0; auto const& mec : mecs) {
        for (auto const& [state, _] : mec) {
            stateToMec[state] = mec_counter;
            mecStates.set(state, true);
        }
        ++mec_counter;
    }

    for (auto const& conjunction : dnf) {
        // get the states of the mdp that (a) are on a MEC, (b) are not already known to be accepting, and (c) don't violate Fins of the conjunction
        storm::storage::BitVector allowed = mecStates & ~acceptingStates;

        if (allowed.empty()) {
            break;  // no more candidates
        }

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
        allMECs += allowedECs.size();
        for (const auto& ec : allowedECs) {
            auto const representativeEcState = ec.begin()->first;
            if (acceptingStates.get(representativeEcState)) {
                // The ec is part of a mec that is already known to be accepting
                continue;
            }

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
                accMECs++;

                // get the MEC containing the current EC
                auto const mec_index = stateToMec[representativeEcState];
                STORM_LOG_ASSERT(mec_index < mecs.size(), "MEC index out of range.");
                auto const& mec = mecs[mec_index];

                // All states of the (global) mec are accepting since we can almost surely reach the inner ec
                for (auto const& [state, _] : mec) {
                    acceptingStates.set(state);
                }

                if (this->isProduceSchedulerSet()) {
                    // save choices for states that weren't assigned to any other MEC yet.
                    this->_schedulerHelper.get().saveProductEcChoices(acceptance, ec, mec, conjunction, product);
                }
            }
        }
    }

    STORM_LOG_INFO("Found " << acceptingStates.getNumberOfSetBits() << " states in " << accMECs << " accepting MECs (considered " << allMECs << " MECs).");

    return acceptingStates;
}

template<typename ValueType, bool Nondeterministic>
storm::storage::BitVector SparseLTLHelper<ValueType, Nondeterministic>::computeAcceptingBCCs(automata::AcceptanceCondition const& acceptance,
                                                                                             storm::storage::SparseMatrix<ValueType> const& transitionMatrix) {
    storm::storage::StronglyConnectedComponentDecomposition<ValueType> bottomSccs(
        transitionMatrix, storage::StronglyConnectedComponentDecompositionOptions().onlyBottomSccs().dropNaiveSccs());
    storm::storage::BitVector acceptingStates(transitionMatrix.getRowGroupCount(), false);

    std::size_t checkedBSCCs = 0, acceptingBSCCs = 0, acceptingBSCCStates = 0;
    for (auto& scc : bottomSccs) {
        checkedBSCCs++;
        if (acceptance.isAccepting(scc)) {
            acceptingBSCCs++;
            for (auto& state : scc) {
                acceptingStates.set(state);
                acceptingBSCCStates++;
            }
        }
    }
    STORM_LOG_INFO("BSCC analysis: " << acceptingBSCCs << " of " << checkedBSCCs << " BSCCs were acceptingStates (" << acceptingBSCCStates
                                     << " states in acceptingStates BSCCs).");
    return acceptingStates;
}

template<typename ValueType, bool Nondeterministic>
std::vector<ValueType> SparseLTLHelper<ValueType, Nondeterministic>::computeDAProductProbabilities(
    Environment const& env, storm::automata::DeterministicAutomaton const& da, std::map<std::string, storm::storage::BitVector>&& apSatSets) {
    auto statesForAP = detail::computeStatesForApSet(da.getAPSet(), std::move(apSatSets));

    storm::storage::BitVector statesOfInterest;
    if (this->hasRelevantStates()) {
        statesOfInterest = this->getRelevantStates();
    } else {
        // Product from all model states
        statesOfInterest = storm::storage::BitVector(this->_transitionMatrix.getRowGroupCount(), true);
    }

    STORM_LOG_INFO("Building " + (Nondeterministic ? std::string("MDP-DA") : std::string("DTMC-DA")) + " product with deterministic automaton, starting from "
                   << statesOfInterest.getNumberOfSetBits() << " model states...");
    transformer::DAProductBuilder productBuilder(da, statesForAP);

    auto product = productBuilder.build<productModelType>(this->_transitionMatrix, statesOfInterest);

    STORM_LOG_INFO("Product " + (Nondeterministic ? std::string("MDP-DA") : std::string("DTMC-DA")) + " has "
                   << product->getProductModel().getNumberOfStates() << " states and " << product->getProductModel().getNumberOfTransitions()
                   << " transitions.");

    // Prepare scheduler
    if (this->isProduceSchedulerSet()) {
        STORM_LOG_THROW(Nondeterministic, storm::exceptions::InvalidOperationException, "Scheduler export only supported for nondeterministic models.");
        this->_schedulerHelper.emplace(product->getProductModel().getNumberOfStates());
    }

    // Compute accepting states
    storm::storage::BitVector acceptingStates;
    if (Nondeterministic) {
        STORM_LOG_INFO("Computing MECs and checking for acceptance...");
        acceptingStates = computeAcceptingECs(*product->getAcceptance(), product->getProductModel().getTransitionMatrix(),
                                              product->getProductModel().getBackwardTransitions(), product);

    } else {
        STORM_LOG_INFO("Computing BSCCs and checking for acceptance...");
        acceptingStates = computeAcceptingBCCs(*product->getAcceptance(), product->getProductModel().getTransitionMatrix());
    }

    if (acceptingStates.empty()) {
        STORM_LOG_INFO("No accepting states, skipping probability computation.");
        if (this->isProduceSchedulerSet()) {
            this->_schedulerHelper.get().setRandom();
        }
        std::vector<ValueType> numericResult(this->_transitionMatrix.getRowGroupCount(), storm::utility::zero<ValueType>());
        return numericResult;
    }

    STORM_LOG_INFO("Computing probabilities for reaching accepting components...");

    storm::storage::BitVector bvTrue(product->getProductModel().getNumberOfStates(), true);
    storm::storage::BitVector soiProduct(product->getStatesOfInterest());

    // Create goal for computeUntilProbabilities, always compute maximizing probabilities
    storm::solver::SolveGoal<ValueType> solveGoalProduct;
    if (this->isValueThresholdSet()) {
        solveGoalProduct = storm::solver::SolveGoal<ValueType>(OptimizationDirection::Maximize, this->getValueThresholdComparisonType(),
                                                               this->getValueThresholdValue(), std::move(soiProduct));
    } else {
        solveGoalProduct = storm::solver::SolveGoal<ValueType>(OptimizationDirection::Maximize);
        solveGoalProduct.setRelevantValues(std::move(soiProduct));
    }

    std::vector<ValueType> prodNumericResult;

    if (Nondeterministic) {
        MDPSparseModelCheckingHelperReturnType<ValueType> prodCheckResult =
            storm::modelchecker::helper::SparseMdpPrctlHelper<ValueType>::computeUntilProbabilities(
                env, std::move(solveGoalProduct), product->getProductModel().getTransitionMatrix(), product->getProductModel().getBackwardTransitions(), bvTrue,
                acceptingStates, this->isQualitativeSet(),
                this->isProduceSchedulerSet()  // Whether to create memoryless scheduler for the Model-DA Product.
            );
        prodNumericResult = std::move(prodCheckResult.values);

        if (this->isProduceSchedulerSet()) {
            this->_schedulerHelper.get().prepareScheduler(da.getNumberOfStates(), acceptingStates, std::move(prodCheckResult.scheduler), productBuilder,
                                                          product, statesOfInterest, this->_transitionMatrix);
        }

    } else {
        prodNumericResult = storm::modelchecker::helper::SparseDtmcPrctlHelper<ValueType>::computeUntilProbabilities(
            env, std::move(solveGoalProduct), product->getProductModel().getTransitionMatrix(), product->getProductModel().getBackwardTransitions(), bvTrue,
            acceptingStates, this->isQualitativeSet());
    }

    std::vector<ValueType> numericResult = product->projectToOriginalModel(this->_transitionMatrix.getRowGroupCount(), prodNumericResult);

    return numericResult;
}

template<typename ValueType, bool Nondeterministic>
std::vector<ValueType> SparseLTLHelper<ValueType, Nondeterministic>::computeLTLProbabilities(Environment const& env, storm::logic::PathFormula const& formula,
                                                                                             std::map<std::string, storm::storage::BitVector>&& apSatSets) {
    std::shared_ptr<storm::logic::Formula const> ltlFormula;
    STORM_LOG_THROW((!Nondeterministic) || this->isOptimizationDirectionSet(), storm::exceptions::InvalidPropertyException,
                    "Formula needs to specify whether minimal or maximal values are to be computed on nondeterministic model.");
    if (Nondeterministic && this->getOptimizationDirection() == OptimizationDirection::Minimize) {
        // negate formula in order to compute 1-Pmax[!formula]
        ltlFormula = std::make_shared<storm::logic::UnaryBooleanPathFormula>(storm::logic::UnaryBooleanOperatorType::Not, formula.asSharedPointer());
        STORM_LOG_INFO("Computing Pmin, proceeding with negated LTL formula.");
    } else {
        ltlFormula = formula.asSharedPointer();
    }

    STORM_LOG_INFO("Resulting LTL path formula: " << ltlFormula->toString());
    STORM_LOG_INFO(" in prefix format: " << ltlFormula->toPrefixString());

    // Convert LTL formula to a deterministic automaton
    STORM_LOG_THROW(env.modelchecker().getLtlAutomatonType() == storm::automata::AutomatonType::DA, storm::exceptions::NotSupportedException,
                    "Only deterministic automata are currently supported for LTL model checking.");
    // For nondeterministic models, the acceptance condition must be in DNF
    auto da = detail::buildAutomatonFromFormula<storm::automata::DeterministicAutomaton>(env, *ltlFormula, Nondeterministic);

    std::vector<ValueType> numericResult = computeDAProductProbabilities(env, *da, std::move(apSatSets));

    if (Nondeterministic && this->getOptimizationDirection() == OptimizationDirection::Minimize) {
        // compute 1-Pmax[!fomula]
        for (auto& value : numericResult) {
            value = storm::utility::one<ValueType>() - value;
        }
    }

    return numericResult;
}
//
// template<typename ValueType, bool Nondeterministic>
// std::tuple<typename SparseLTLHelper<ValueType, Nondeterministic>::productModelType, std::vector<storm::automata::AcceptanceCondition::ptr>,
//           std::vector<uint64_t>>
// SparseLTLHelper<ValueType, Nondeterministic>::buildFromFormulas(productModelType const& model,
//                                                                std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas,
//                                                                Environment const& env) {
//    std::vector<storm::automata::AcceptanceCondition::ptr> acceptanceConditions(formulas.size());
//    productModelType productModel(model.getTransitionMatrix(), model.getStateLabeling());
//    std::vector<storm::storage::BitVector> modelStateToProductStates(model.getNumberOfStates());
//    for (uint64_t i = 0; i < model.getNumberOfStates(); i++) {
//        auto iAsVector = storage::BitVector(model.getNumberOfStates());
//        iAsVector.set(i);
//        modelStateToProductStates[i] = iAsVector;
//    }
//
//    for (int i = 0; i < formulas.size(); i++) {
//        // compute product
//        auto const& pathformula = formulas[i]->asOperatorFormula().getSubformula().asPathFormula();
//        auto product = buildFromFormula(productModel, pathformula, env);
//        acceptanceConditions[i] = product->getAcceptance();
//
//        // lift state labeling
//        models::sparse::StateLabeling productLabeling(product->getProductModel().getNumberOfStates());
//        for (auto const& label : productModel.getStateLabeling().getLabels()) {
//            productLabeling.addLabel(label);
//            if (label == "init")
//                continue;
//
//            auto statesModel = productModel.getStateLabeling().getStates(label);
//            auto statesProduct = product->liftFromModel(statesModel);
//            productLabeling.setStates(label, statesProduct);
//        }
//
//        productLabeling.setStates("init", product->getProductModel().getStateLabeling().getStates("soi"));
//        productModel = productModelType(product->getProductModel().getTransitionMatrix(), productLabeling);
//
//        // lift state mapping
//        for (int j = 0; j < model.getNumberOfStates(); j++) {
//            modelStateToProductStates[j] = product->liftFromModel(modelStateToProductStates[j]);
//        }
//
//        for (int j = 0; j < i; j++) {
//            // lift product acceptance condition
//            acceptanceConditions[j] = acceptanceConditions[j]->lift(product->getProductModel().getNumberOfStates(),
//                                                                    [&product](std::size_t prodState) { return product->getModelState(prodState); });
//        }
//    }
//
//    std::vector<uint64_t> indexToModelState(productModel.getNumberOfStates());
//    for (uint64_t i = 0; i < model.getNumberOfStates(); i++) {
//        for (auto const j : modelStateToProductStates[i]) {
//            indexToModelState[j] = i;
//        }
//    }
//
//    return std::make_tuple(productModel, acceptanceConditions, indexToModelState);
//}

template<typename SparseModelType>
ToRabinReturnType<SparseModelType> LTLToRabinObjectives(Environment const& env, SparseModelType const& originalModel,
                                                        std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
                                                        std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& formulaChecker) {
    using ValueType = typename SparseModelType::ValueType;
    bool constexpr Nondeterministic =
        std::is_same_v<SparseModelType, models::sparse::Mdp<ValueType>> || std::is_same_v<SparseModelType, models::sparse::MarkovAutomaton<ValueType>>;
    STORM_LOG_ASSERT(Nondeterministic == originalModel.isNondeterministicModel(), "Unexpected model type.");

    ToRabinReturnType<SparseModelType> result;
    auto productToOriginalMap = storm::utility::vector::buildVectorForRange<uint64_t>(0, originalModel.getNumberOfStates());  // start with identity mapping

    for (auto const& f : ltlFormulas) {
        STORM_LOG_THROW(f->isPathFormula(), storm::exceptions::InvalidPropertyException, "Expected path formula, got " << *f << ".");
        // Replace state-subformulae by atomic propositions using PCTL*-style model checking
        storm::logic::ExtractMaximalStateFormulasVisitor::ApToFormulaMap extracted;
        auto ltlFormula = storm::logic::ExtractMaximalStateFormulasVisitor::extract(f->asPathFormula(), extracted);
        STORM_LOG_ASSERT(ltlFormula->isPathFormula(), "Unexpected formula type.");
        auto apSatSetsOrigModel = SparseLTLHelper<ValueType, Nondeterministic>::computeApSets(extracted, formulaChecker);
        auto const& modelRef = result.model == nullptr ? originalModel : *result.model;

        typename transformer::DAProduct<SparseModelType>::ptr product;
        using AutomatonType = storm::automata::AutomatonType;
        if (env.modelchecker().getLtlAutomatonType() == AutomatonType::LDBA) {
            auto ldba = detail::buildAutomatonFromFormula<storm::automata::LimitDeterministicAutomaton>(env, *ltlFormula, true);  // always force DNF
            auto statesForAP = detail::computeStatesForApSet(ldba->getAPSet(), std::move(apSatSetsOrigModel), productToOriginalMap);
            assert(false);  // todo
            // transformer::LDBAProductBuilder productBuilder(*automaton, statesForAP);
            // product = productBuilder.build(modelRef);
        } else {
            STORM_LOG_ASSERT(env.modelchecker().getLtlAutomatonType() == AutomatonType::DA, "Unsupported automaton type.");
            auto da = detail::buildAutomatonFromFormula<storm::automata::DeterministicAutomaton>(env, *ltlFormula, true);  // always force DNF
            auto statesForAP = detail::computeStatesForApSet(da->getAPSet(), std::move(apSatSetsOrigModel), productToOriginalMap);
            assert(false);  // todo
            // transformer::DAProductBuilder productBuilder(*automaton, statesForAP);
            // product = productBuilder.build(modelRef);
        }
        productToOriginalMap = product->liftFromModel(productToOriginalMap);
        result.model = product->getProductModel().template as<SparseModelType>();
        result.rabinObjectives.push_back(
            detail::extractRabinObjectiveFromProductAcceptance(std::move(*product->getAcceptance()), result.model->getStateLabeling()));
    }
    return result;
}

template class SparseLTLHelper<double, false>;
template class SparseLTLHelper<double, true>;

template ToRabinReturnType<storm::models::sparse::Mdp<double>> LTLToRabinObjectives<storm::models::sparse::Mdp<double>>(
    Environment const& env, storm::models::sparse::Mdp<double> const& originalModel,
    std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
    std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& formulaChecker);
template ToRabinReturnType<storm::models::sparse::MarkovAutomaton<double>> LTLToRabinObjectives<storm::models::sparse::MarkovAutomaton<double>>(
    Environment const& env, storm::models::sparse::MarkovAutomaton<double> const& originalModel,
    std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
    std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& formulaChecker);

#ifdef STORM_HAVE_CARL
template class SparseLTLHelper<storm::RationalNumber, false>;
template class SparseLTLHelper<storm::RationalNumber, true>;
template class SparseLTLHelper<storm::RationalFunction, false>;

template ToRabinReturnType<storm::models::sparse::Mdp<storm::RationalNumber>> LTLToRabinObjectives<storm::models::sparse::Mdp<storm::RationalNumber>>(
    Environment const& env, storm::models::sparse::Mdp<storm::RationalNumber> const& originalModel,
    std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
    std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& formulaChecker);
template ToRabinReturnType<storm::models::sparse::MarkovAutomaton<storm::RationalNumber>>
LTLToRabinObjectives<storm::models::sparse::MarkovAutomaton<storm::RationalNumber>>(
    Environment const& env, storm::models::sparse::MarkovAutomaton<storm::RationalNumber> const& originalModel,
    std::vector<std::shared_ptr<storm::logic::Formula const>> const& ltlFormulas,
    std::function<storm::storage::BitVector(storm::logic::Formula const&)> const& formulaChecker);

#endif

}  // namespace helper
}  // namespace modelchecker
}  // namespace storm

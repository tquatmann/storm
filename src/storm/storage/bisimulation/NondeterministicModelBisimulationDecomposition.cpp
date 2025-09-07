#include "storm/storage/bisimulation/NondeterministicModelBisimulationDecomposition.h"

#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"

#include "storm/utility/graph.h"

#include "storm/exceptions/IllegalFunctionCallException.h"
#include "storm/utility/macros.h"

#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace storage {

using namespace bisimulation;

template<typename ModelType>
NondeterministicModelBisimulationDecomposition<ModelType>::NondeterministicModelBisimulationDecomposition(
    ModelType const& model,
    typename BisimulationDecomposition<ModelType>::Options const& options)
    : BisimulationDecomposition<ModelType>(model, model.getTransitionMatrix().transpose(false),
                                                                                                          options),
      choiceToStateMapping(model.getNumberOfChoices()),
      quotientDistributions(model.getNumberOfChoices()),
      orderedQuotientDistributions(model.getNumberOfChoices()) {
    STORM_LOG_THROW(options.getType() == BisimulationType::Strong, storm::exceptions::IllegalFunctionCallException,
                    "Weak bisimulation is currently not supported for nondeterministic models.");
}

template<typename ModelType>
std::pair<storm::storage::BitVector, storm::storage::BitVector> NondeterministicModelBisimulationDecomposition<ModelType>::getStatesWithProbability01() {
    STORM_LOG_THROW(this->options.isOptimizationDirectionSet(), storm::exceptions::IllegalFunctionCallException,
                    "Can only compute states with probability 0/1 with an optimization direction (min/max).");
    if (this->options.getOptimizationDirection() == OptimizationDirection::Minimize) {
        return storm::utility::graph::performProb01Min(this->model.getTransitionMatrix(), this->model.getTransitionMatrix().getRowGroupIndices(),
                                                       this->model.getBackwardTransitions(), this->options.phiStates.get(), this->options.psiStates.get());
    } else {
        return storm::utility::graph::performProb01Max(this->model.getTransitionMatrix(), this->model.getTransitionMatrix().getRowGroupIndices(),
                                                       this->model.getBackwardTransitions(), this->options.phiStates.get(), this->options.psiStates.get());
    }
}

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::initialize() {
    this->createChoiceToStateMapping();
    this->initializeQuotientDistributions();
}

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::createChoiceToStateMapping() {
    std::vector<uint_fast64_t> nondeterministicChoiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    for (storm::storage::sparse::state_type state = 0; state < this->model.getNumberOfStates(); ++state) {
        for (uint_fast64_t choice = nondeterministicChoiceIndices[state]; choice < nondeterministicChoiceIndices[state + 1]; ++choice) {
            choiceToStateMapping[choice] = state;
        }
    }
}

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::initializeQuotientDistributions() {
    // TODO: implement
}

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::updateOrderedQuotientDistributions(storm::storage::sparse::state_type state) {
    std::vector<uint_fast64_t> nondeterministicChoiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    std::sort(this->orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state],
              this->orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state + 1],
              [this](storm::storage::Distribution<ValueType> const* dist1, storm::storage::Distribution<ValueType> const* dist2) {
                  return dist1->less(*dist2, this->comparator);
              });
}

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::buildQuotient() {
    // TODO: implement
}

template<typename ModelType>
bool NondeterministicModelBisimulationDecomposition<ModelType>::checkQuotientDistributions() const {
    // TODO: implement
}

template<typename ModelType>
bool NondeterministicModelBisimulationDecomposition<ModelType>::printDistributions(uint_fast64_t state) const {
    std::vector<uint_fast64_t> nondeterministicChoiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();
    for (auto choice = nondeterministicChoiceIndices[state]; choice < nondeterministicChoiceIndices[state + 1]; ++choice) {
        std::cout << quotientDistributions[choice] << '\n';
    }
    return true;
}

// template<typename ModelType>
// bool NondeterministicModelBisimulationDecomposition<ModelType>::splitBlockAccordingToCurrentQuotientDistributions() {
//     // TODO: implement
// }

template<typename ModelType>
bool NondeterministicModelBisimulationDecomposition<ModelType>::quotientDistributionsLess(storm::storage::sparse::state_type state1,
                                                                                          storm::storage::sparse::state_type state2) const {
    STORM_LOG_TRACE("Comparing the quotient distributions of state " << state1 << " and " << state2 << ".");
    std::vector<uint_fast64_t> nondeterministicChoiceIndices = this->model.getTransitionMatrix().getRowGroupIndices();

    auto firstIt = orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state1];
    auto firstIte = orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state1 + 1];
    auto secondIt = orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state2];
    auto secondIte = orderedQuotientDistributions.begin() + nondeterministicChoiceIndices[state2 + 1];

    for (; firstIt != firstIte && secondIt != secondIte; ++firstIt, ++secondIt) {
        // If the current distributions are in a less-than relationship, we can return a result.
        if ((*firstIt)->less(**secondIt, this->comparator)) {
            return true;
        } else if ((*secondIt)->less(**firstIt, this->comparator)) {
            return false;
        }

        // If the distributions matched, we need to advance both distribution iterators to the next distribution
        // that is larger.
        while (firstIt != firstIte && std::next(firstIt) != firstIte && !(*firstIt)->less(**std::next(firstIt), this->comparator)) {
            ++firstIt;
        }
        while (secondIt != secondIte && std::next(secondIt) != secondIte && !(*secondIt)->less(**std::next(secondIt), this->comparator)) {
            ++secondIt;
        }
    }

    if (firstIt == firstIte && secondIt != secondIte) {
        return true;
    }
    return false;
}

// template<typename ModelType>
// void NondeterministicModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(
// // TODO: implement
// }

template<typename ModelType>
void NondeterministicModelBisimulationDecomposition<ModelType>::refinePartitionBasedOnSplitter(std::span<uint64_t const> splitterBlock, std::deque<typename bisimulation::Partition::Block>& splitterQueue,
                                                                                               bisimulation::Partition::BlockSet& enqueuedSplitterBlocks) {

}

template class NondeterministicModelBisimulationDecomposition<storm::models::sparse::Mdp<double>>;

#ifdef STORM_HAVE_CARL
template class NondeterministicModelBisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalNumber>>;
template class NondeterministicModelBisimulationDecomposition<storm::models::sparse::Mdp<storm::RationalFunction>>;

template class storm::storage::NondeterministicModelBisimulationDecomposition<
    storm::models::sparse::Mdp<carl::Interval<double>>>;
#endif
}  // namespace storage
}  // namespace storm

#include "storm/modelchecker/prctl/helper/BaierUpperRewardBoundsComputer.h"

#include "storm/adapters/RationalNumberAdapter.h"

#include "storm/storage/BitVector.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/StronglyConnectedComponentDecomposition.h"

#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

namespace storm {
namespace modelchecker {
namespace helper {

template<typename ValueType>
BaierUpperRewardBoundsComputer<ValueType>::BaierUpperRewardBoundsComputer(storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                                          std::vector<ValueType> const& rewards,
                                                                          std::vector<ValueType> const& oneStepTargetProbabilities,
                                                                          std::function<uint64_t(uint64_t)> const& stateToScc)
    : _transitionMatrix(transitionMatrix), _stateToScc(stateToScc), _rewards(rewards), _oneStepTargetProbabilities(oneStepTargetProbabilities) {
    // Intentionally left empty.
}

template<typename ValueType>
std::vector<ValueType> BaierUpperRewardBoundsComputer<ValueType>::computeUpperBoundOnExpectedVisitingTimes(
    storm::storage::SparseMatrix<ValueType> const& transitionMatrix, std::vector<ValueType> const& oneStepTargetProbabilities) {
    std::vector<uint64_t> stateToScc(transitionMatrix.getRowGroupCount());
    {
        // Create an SCC decomposition of the system.
        storm::storage::StronglyConnectedComponentDecomposition<ValueType> sccDecomposition(transitionMatrix);

        uint64_t sccIndex = 0;
        for (auto const& block : sccDecomposition) {
            for (auto const& state : block) {
                stateToScc[state] = sccIndex;
            }
            ++sccIndex;
        }
    }
    return computeUpperBoundOnExpectedVisitingTimes(transitionMatrix, oneStepTargetProbabilities, [&stateToScc](uint64_t s) { return stateToScc[s]; });
}

template<typename ValueType>
std::vector<ValueType> BaierUpperRewardBoundsComputer<ValueType>::computeUpperBoundOnExpectedVisitingTimes(
    storm::storage::SparseMatrix<ValueType> const& transitionMatrix, std::vector<ValueType> const& oneStepTargetProbabilities,
    std::function<uint64_t(uint64_t)> const& stateToScc) {
    assert(transitionMatrix.getRowCount() == oneStepTargetProbabilities.size());
    auto const numStates = transitionMatrix.getRowGroupCount();

    // A choice is valid iff it goes to non-remaining states with non-zero probability.
    // Initially, mark all choices as valid that have non-zero probability to go to the target states *or* to a different Scc.
    auto validChoices = storm::utility::vector::filterGreaterZero(oneStepTargetProbabilities);
    for (uint64_t state = 0; state < numStates; ++state) {
        auto const scc = stateToScc(state);
        for (auto rowIndex = transitionMatrix.getRowGroupIndices()[state], rowEnd = transitionMatrix.getRowGroupIndices()[state + 1]; rowIndex < rowEnd;
             ++rowIndex) {
            auto const row = transitionMatrix.getRow(rowIndex);
            if (std::any_of(row.begin(), row.end(), [&stateToScc, &scc](auto const& entry) { return scc != stateToScc(entry.getColumn()); })) {
                validChoices.set(rowIndex, true);
            }
        }
    }

    // Vector that holds the result.
    std::vector<ValueType> result(numStates, storm::utility::one<ValueType>());

    // The states that we still need to assign a value.
    storm::storage::BitVector processedStates(numStates, false);

    auto backwardTransitions = transitionMatrix.transpose(true);

    // "Over-approximates" set S_i from the paper
    storm::storage::BitVector currentStates(numStates, true);
    // "Over-approximates" set S_{i+1}
    storm::storage::BitVector nextStates(numStates, false);
    // We usually get the best performance by processing states in inverted order.
    // This is because in the sparse engine the states are explored with a BFS/DFS
    // To achive processing in a backwards order, the current and next states are stored with their "inverted" indices.
    auto const largestIndex = numStates - 1;
    while (true) {
        for (auto stateInv : currentStates) {
            auto state = largestIndex - stateInv;
            if (processedStates.get(state)) {
                continue;  // do not process already processed states.
            }
            auto const& startRow = transitionMatrix.getRowGroupIndices()[state];
            auto const& endRow = transitionMatrix.getRowGroupIndices()[state + 1];
            // update the valid choices at the current states
            bool canProcessState = true;
            // Check if state belongs to the set "S_i" from the paper
            for (auto rowIndex = validChoices.getNextUnsetIndex(startRow); rowIndex < endRow; rowIndex = validChoices.getNextUnsetIndex(rowIndex + 1)) {
                auto row = transitionMatrix.getRow(rowIndex);
                if (std::any_of(row.begin(), row.end(), [&processedStates](auto const& entry) { return processedStates.get(entry.getColumn()); })) {
                    validChoices.set(rowIndex, true);
                } else {
                    canProcessState = false;
                    break;
                }
            }
            if (canProcessState) {
                // The state is in the set S_i, i.e. we can compute a value for it.
                auto const scc = stateToScc(state);
                auto& stateValue = result[state];
                for (auto rowIndex = startRow; rowIndex < endRow; ++rowIndex) {
                    ValueType rowValue = oneStepTargetProbabilities[rowIndex];
                    for (auto const& entry : transitionMatrix.getRow(rowIndex)) {
                        if (auto successorState = entry.getColumn(); scc != stateToScc(successorState)) {
                            rowValue += entry.getValue();  // * 1
                        } else if (processedStates.get(successorState)) {
                            rowValue += entry.getValue() * result[successorState];
                        }
                    }
                    stateValue = std::min(stateValue, rowValue);
                }
                processedStates.set(state);

                // Iterate through predecessors that might become valid in the next run
                for (auto const& predEntry : backwardTransitions.getRow(state)) {
                    nextStates.set(largestIndex - predEntry.getColumn());  // Insert inverted index
                }
            }
        }
        if (nextStates.empty()) {
            break;
        }
        currentStates.clear();
        std::swap(nextStates, currentStates);
    }

    // Transform the d_t to an upper bound for zeta(t)
    storm::utility::vector::applyPointwise(result, result, [](ValueType const& r) -> ValueType { return storm::utility::one<ValueType>() / r; });
    return result;
}

template<typename ValueType>
ValueType BaierUpperRewardBoundsComputer<ValueType>::computeUpperBound() {
    STORM_LOG_TRACE("Computing upper reward bounds using variant-2 of Baier et al.");
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    auto expVisits = _stateToScc ? computeUpperBoundOnExpectedVisitingTimes(_transitionMatrix, _oneStepTargetProbabilities, _stateToScc)
                                 : computeUpperBoundOnExpectedVisitingTimes(_transitionMatrix, _oneStepTargetProbabilities);

    ValueType upperBound = storm::utility::zero<ValueType>();
    for (uint64_t state = 0; state < expVisits.size(); ++state) {
        ValueType maxReward = storm::utility::zero<ValueType>();
        // By starting the maxReward with zero, negative rewards are essentially ignored which
        // is necessary to provide a valid upper bound
        for (auto row = _transitionMatrix.getRowGroupIndices()[state], endRow = _transitionMatrix.getRowGroupIndices()[state + 1]; row < endRow; ++row) {
            maxReward = std::max(maxReward, _rewards[row]);
        }
        upperBound += expVisits[state] * maxReward;
    }

    STORM_LOG_TRACE("Baier algorithm for reward bound computation (variant 2) computed bound " << upperBound << ".");
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    STORM_LOG_TRACE("Computed upper bounds on rewards in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms.");
    return upperBound;
}

template class BaierUpperRewardBoundsComputer<double>;

#ifdef STORM_HAVE_CARL
template class BaierUpperRewardBoundsComputer<storm::RationalNumber>;
#endif
}  // namespace helper
}  // namespace modelchecker
}  // namespace storm

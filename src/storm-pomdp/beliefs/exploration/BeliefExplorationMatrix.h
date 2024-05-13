#pragma once

#include <tuple>
#include <vector>

#include "storm-pomdp/beliefs/utility/types.h"

namespace storm::pomdp::beliefs {

template<typename ValueType, typename...>
struct BeliefExplorationTransition {
    ValueType probability;
    storm::pomdp::beliefs::BeliefId targetBelief;
};

template<typename ValueType, typename ExtraTransitionData>
struct BeliefExplorationTransition<ValueType, ExtraTransitionData> {
    ValueType probability;
    storm::pomdp::beliefs::BeliefId targetBelief;
    ExtraTransitionData data;
};

template<typename ValueType, typename FirstExtraTransitionData, typename... OtherExtraTransitionData>
struct BeliefExplorationTransition<ValueType, FirstExtraTransitionData, OtherExtraTransitionData...> {
    ValueType probability;
    storm::pomdp::beliefs::BeliefId targetBelief;
    std::tuple<FirstExtraTransitionData, OtherExtraTransitionData...> data;
};

template<typename ValueType, typename... ExtraTransitionData>
class BeliefExplorationMatrix {
   public:
    /*!
     * Initializes a new (empty) belief exploration matrix.
     */
    BeliefExplorationMatrix();

    /*!
     * While building the matrix, ends the current row in the matrix.
     */
    void endCurrentRow();

    /*!
     * While building the matrix, ends the current row group in the matrix.
     * @note This function should be called after endCurrentRow() has been called.
     */
    void endCurrentRowGroup();

    /*!
     * @return the current number of rows in the matrix
     */
    std::size_t rows() const;

    /*!
     * @return the current number of row groups in the matrix
     */
    std::size_t groups() const;

    /*!
     * Stores the successor
     */
    std::vector<BeliefExplorationTransition<ValueType, ExtraTransitionData...>> transitions;
    std::vector<uint64_t> rowIndications;
    std::vector<uint64_t> rowGroupIndices;
};

}  // namespace storm::pomdp::beliefs
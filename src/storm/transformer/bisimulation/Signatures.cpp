#include "storm/transformer/bisimulation/Signatures.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

template<typename StateSignature>
Signatures<StateSignature>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                       storm::bisimulation::Partition const& partition)
    : model(model), partition(partition), choiceClasses(choiceClasses), stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType>
bool ExactStateSignature<ValueType>::ChoiceSignature::operator<(ChoiceSignature const& other) const {
    if (choiceClass != other.choiceClass) {
        return choiceClass < other.choiceClass;
    }
    if (distr.size() != other.distr.size()) {
        return distr.size() < other.distr.size();
    }
    auto it2 = other.distr.begin();
    for (auto const& entry : distr) {
        if (entry.first.size() != it2->first.size()) {
            return entry.first.size() < it2->first.size();
        }
        if (entry.first.data() != it2->first.data()) {
            return entry.first.data() < it2->first.data();
        }
        if (entry.second != it2->second) {
            return entry.second < it2->second;
        }
        ++it2;
    }
    return false;
}

template<typename StateSignature>
void Signatures<StateSignature>::updateStateSignature(uint64_t const stateIndex) {
    static_assert(std::is_same_v<StateSignature, ExactStateSignature<typename StateSignature::ValueType>>,
                  "Only exact state signatures are supported for now.");

    StateSignature& signature = stateSignatureCache[stateIndex];
    signature.choices.clear();

    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        typename StateSignature::ChoiceSignature choice{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0, .distr = {}};
        for (auto const& entry : model.getTransitionMatrix().getRow(choiceIndex)) {
            if (!storm::utility::isZero(entry.getValue())) {
                auto const emplace_res = choice.distr.emplace(partition.getBlockOfElement(entry.getColumn()), entry.getValue());
                if (!emplace_res.second) {
                    emplace_res.first->second += entry.getValue();
                }
            }
        }
        signature.choices.insert(std::move(choice));
    }
}

template<typename StateSignature>
SplitOrder<StateSignature> Signatures<StateSignature>::getExactSplitOrder() const
    requires IsExact
{
    return SplitOrder<StateSignature>{*this};
}

template<typename StateSignature>
bool SplitOrder<StateSignature>::operator()(uint64_t const state1, uint64_t const state2) const {
    auto const& sig1 = signatures.stateSignatureCache[state1];
    auto const& sig2 = signatures.stateSignatureCache[state2];
    return sig1.choices < sig2.choices;
}

template struct ExactStateSignature<double>;
template struct ExactStateSignature<storm::RationalNumber>;
template struct ExactStateSignature<storm::RationalFunction>;
template struct ExactStateSignature<storm::Interval>;
template struct ExactStateSignature<storm::RationalInterval>;

template class Signatures<ExactStateSignature<double>>;
template class Signatures<ExactStateSignature<storm::RationalNumber>>;
template class Signatures<ExactStateSignature<storm::RationalFunction>>;
template class Signatures<ExactStateSignature<storm::Interval>>;
template class Signatures<ExactStateSignature<storm::RationalInterval>>;

template struct SplitOrder<ExactStateSignature<double>>;
template struct SplitOrder<ExactStateSignature<storm::RationalNumber>>;
template struct SplitOrder<ExactStateSignature<storm::RationalFunction>>;
template struct SplitOrder<ExactStateSignature<storm::Interval>>;
template struct SplitOrder<ExactStateSignature<storm::RationalInterval>>;

}  // namespace storm::bisimulation

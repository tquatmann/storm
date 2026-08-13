#include "storm/transformer/bisimulation/Signature.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/constants.h"

namespace storm::bisimulation {

template<typename ValueType>
bool Signatures<ValueType>::ChoiceSignature::operator<(ChoiceSignature const& other) const {
    if (choiceClass != other.choiceClass) {
        return choiceClass < other.choiceClass;
    }
    // std::map's relational operators are unavailable here since Block (a std::span) has no operator</operator==,
    // so we compare manually using BlockCompare (both maps are already ordered by it).
    Partition::BlockCompare const blockLess;
    auto it1 = distr.begin();
    auto it2 = other.distr.begin();
    for (; it1 != distr.end() && it2 != other.distr.end(); ++it1, ++it2) {
        if (blockLess(it1->first, it2->first)) {
            return true;
        }
        if (blockLess(it2->first, it1->first)) {
            return false;
        }
        if (it1->second != it2->second) {
            return it1->second < it2->second;
        }
    }
    return it1 == distr.end() && it2 != other.distr.end();
}

template<typename ValueType>
bool Signatures<ValueType>::ChoiceSignature::operator==(ChoiceSignature const& other) const {
    if (choiceClass != other.choiceClass || distr.size() != other.distr.size()) {
        return false;
    }
    Partition::BlockCompare const blockLess;
    auto it2 = other.distr.begin();
    for (auto it1 = distr.begin(); it1 != distr.end(); ++it1, ++it2) {
        if (blockLess(it1->first, it2->first) || blockLess(it2->first, it1->first) || it1->second != it2->second) {
            return false;
        }
    }
    return true;
}

template<typename ValueType>
bool Signatures<ValueType>::StateSignature::operator<(StateSignature const& other) const {
    return choices < other.choices;
}

template<typename ValueType>
Signatures<ValueType>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                  storm::bisimulation::Partition const& partition)
    : model(model), partition(partition), choiceClasses(choiceClasses), stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType>
auto Signatures<ValueType>::computeChoiceSignature(uint64_t const choiceIndex) const -> ChoiceSignature {
    Partition::OrderedBlockMap<ValueType> distr;
    for (auto const& entry : model.getTransitionMatrix().getRow(choiceIndex)) {
        if (!storm::utility::isZero(entry.getValue())) {
            auto const emplace_res = distr.emplace(partition.getBlockOfElement(entry.getColumn()), entry.getValue());
            if (!emplace_res.second) {
                emplace_res.first->second += entry.getValue();
            }
        }
    }
    return ChoiceSignature{choiceClasses ? (*choiceClasses)[choiceIndex] : 0, std::move(distr)};
}

template<typename ValueType>
void Signatures<ValueType>::updateStateSignature(uint64_t const stateIndex) {
    StateSignature signature;
    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        signature.choices.insert(computeChoiceSignature(choiceIndex));
    }
    stateSignatureCache[stateIndex] = std::move(signature);
}

template<typename ValueType>
auto Signatures<ValueType>::getStateSignature(uint64_t const stateIndex) const -> StateSignature const& {
    return stateSignatureCache[stateIndex];
}

template class Signatures<double>;
template class Signatures<storm::RationalNumber>;
template class Signatures<storm::RationalFunction>;
template class Signatures<storm::Interval>;
template class Signatures<storm::RationalInterval>;

}  // namespace storm::bisimulation

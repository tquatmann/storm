#include "storm/transformer/bisimulation/Signatures.h"

#include <algorithm>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

template<typename ValueType>
ChoiceSignatureCache<ValueType>::ChoiceSignatureCache(storm::bisimulation::Partition const& partition)
    : values(partition.getNumberOfElements(), storm::utility::zero<ValueType>()) {}

template<typename ValueType>
void ChoiceSignatureCache<ValueType>::addValue(Partition::Block const& b, ValueType const& value) {
    auto& current = values[b.front()];
    if (storm::utility::isZero(current)) {
        support.push_back(b);
    }
    current += value;
}

template<typename ValueType>
std::span<std::pair<Partition::Block, ValueType>> ChoiceSignatureCache<ValueType>::extract(std::pair<Partition::Block, ValueType>* dest) {
    std::sort(support.begin(), support.end(), Partition::BlockCompare());
    auto* const begin = dest;
    for (auto const& b : support) {
        auto& value = values[b.front()];
        *dest = std::make_pair(b, std::move(value));
        ++dest;
        value = storm::utility::zero<ValueType>();
    }
    std::span<std::pair<Partition::Block, ValueType>> const result(begin, support.size());
    support.clear();
    return result;
}

template<typename ValueType, SignatureMode Mode>
std::strong_ordering StateSignature<ValueType, Mode>::ChoiceSignature::compare(ChoiceSignature const& other, [[maybe_unused]] ValueType const tolerance) const {
    // Let c_1, c_2, c_3 be choices and let ≈ be the used equality for the mode (operator== for Exact, approxEqual for Approximate).
    // It must hold that c_1 < c_2 < c_3 and c_1 ≈ c_3 implies c_1 ≈ c_2 ≈ c_3.
    // Since in all cases ≈ only holds if choiceClass and the support of distr coincide, we check those first
    if (auto const cmp = choiceClass <=> other.choiceClass; cmp != 0) {
        return cmp;
    }
    if (auto const cmp = distr.size() <=> other.distr.size(); cmp != 0) {
        return cmp;
    }
    auto it2 = other.distr.begin();
    for (auto const& entry : distr) {
        if (auto const cmp = entry.first.size() <=> it2->first.size(); cmp != 0) {
            return cmp;
        }
        if (auto const cmp = entry.first.data() <=> it2->first.data(); cmp != 0) {
            return cmp;
        }
        ++it2;
    }
    // At this point, only the distribution values could potentially differ.
    // In approximative mode, it is important that these are checked last and that we do it such that approximately equivalent distributions will be close in a
    // ordered sequence.
    // ValueType (e.g. RationalFunction, Interval) does not necessarily support <=>, so we fall back to </!= here.
    std::strong_ordering cmp = std::strong_ordering::equal;
    it2 = other.distr.begin();
    for (auto const& entry : distr) {
        auto compareValues = [&entry, &it2]() { return entry.second < it2->second ? std::strong_ordering::less : std::strong_ordering::greater; };
        if constexpr (Mode == SignatureMode::Exact) {
            if (entry.second != it2->second) {
                return compareValues();
            }
        } else {
            if (storm::utility::abs<ValueType>(entry.second - it2->second) > tolerance) {
                return compareValues();
            } else if (cmp == std::strong_ordering::equal && entry.second != it2->second) {
                cmp = compareValues();
            }
        }
        ++it2;
    }
    // Reaching this point means that all values are equal (in approximative mode: up to the given tolerance).
    // In exact mode, we always return equal. In approximative mode, we return the comparison of the first unequal (without tolerance) pair of values (if any).
    return cmp;
}

template<typename ValueType, SignatureMode Mode>
bool StateSignature<ValueType, Mode>::ChoiceSignature::operator==(ChoiceSignature const& other) const
    requires(Mode == SignatureMode::Exact)
{
    return compare(other, storm::utility::zero<ValueType>()) == std::strong_ordering::equal;
}

template<typename ValueType, SignatureMode Mode>
bool StateSignature<ValueType, Mode>::ChoiceSignature::approxEqual(ChoiceSignature const& other, ValueType const tolerance) const
    requires(Mode == SignatureMode::Approximative)
{
    if (choiceClass != other.choiceClass || distr.size() != other.distr.size()) {
        return false;
    }
    auto it2 = other.distr.begin();
    for (auto const& entry : distr) {
        if (entry.first.data() != it2->first.data() || entry.first.size() != it2->first.size()) {
            return false;
        }
        if (storm::utility::abs<ValueType>(entry.second - it2->second) > tolerance) {
            return false;
        }
        ++it2;
    }
    return true;
}

template<typename ValueType, SignatureMode Mode>
Signatures<ValueType, Mode>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                        storm::bisimulation::Partition const& partition)
    requires(Mode == SignatureMode::Exact)
    : model(model),
      partition(partition),
      choiceClasses(choiceClasses),
      choiceSignatureCache(partition),
      choiceDistributionStorage(model.getTransitionMatrix().getEntryCount()),
      halfTolerance(storm::utility::zero<ValueType>()),
      stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
Signatures<ValueType, Mode>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                        storm::bisimulation::Partition const& partition, ValueType const& tolerance)
    requires(Mode == SignatureMode::Approximative)
    : model(model),
      partition(partition),
      choiceClasses(choiceClasses),
      choiceSignatureCache(partition),
      choiceDistributionStorage(model.getTransitionMatrix().getEntryCount()),
      halfTolerance(tolerance / storm::utility::convertNumber<ValueType, uint64_t>(2)),
      stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getChoiceSignature(uint64_t const choiceIndex) const -> ChoiceSignature {
    auto const& matrix = model.getTransitionMatrix();
    // Use the choiceSignatureCache to compute the distribution over successor blocks
    for (auto const& entry : matrix.getRow(choiceIndex)) {
        if (!storm::utility::isZero(entry.getValue())) {
            choiceSignatureCache.addValue(partition.getBlockOfElement(entry.getColumn()), entry.getValue());
        }
    }
    // Extract the resulting distribution from the cache. We put it into
    uint64_t const offset = std::distance(matrix.begin(), matrix.begin(choiceIndex));
    return ChoiceSignature{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0,
                           .distr = choiceSignatureCache.extract(choiceDistributionStorage.data() + offset)};
}

template<typename ValueType, SignatureMode Mode>
void Signatures<ValueType, Mode>::updateStateSignature(uint64_t const stateIndex) {
    auto& choices = stateSignatureCache[stateIndex].choices;
    choices.clear();
    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        auto choiceSignature = getChoiceSignature(choiceIndex);
        auto [it, found] = stateSignatureCache[stateIndex].find(choiceSignature, halfTolerance);
        if (!found) {
            choices.insert(it, std::move(choiceSignature));
        }
    }
}

template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getSplitOrder() const -> SplitOrder {
    return SplitOrder{.signatures = stateSignatureCache, .tolerance = halfTolerance};
}

template<typename ValueType, SignatureMode Mode>
bool Signatures<ValueType, Mode>::SplitOrder::operator()(uint64_t const state1, uint64_t const state2) const {
    auto const& sig1 = signatures[state1];
    auto const& sig2 = signatures[state2];
    // Do not compare the originating choice index (i.e. only compare the set of keys of the respective choices maps)
    if (auto const cmp = sig1.choices.size() <=> sig2.choices.size(); cmp != 0) {
        return cmp < 0;
    }
    auto it2 = sig2.choices.begin();
    for (auto const& choice : sig1.choices) {
        if (auto const cmp = choice.compare(*it2, tolerance); cmp != 0) {
            return cmp < 0;
        }
        ++it2;
    }
    return false;
}
template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getSplitCondition() const -> SplitCondition {
    return SplitCondition{.signatures = stateSignatureCache, .tolerance = halfTolerance};
}

template<typename ValueType, SignatureMode Mode>
bool Signatures<ValueType, Mode>::SplitCondition::operator()(uint64_t const state1, uint64_t const state2) const {
    auto const& sig1 = signatures[state1];
    auto const& sig2 = signatures[state2];
    if (sig1.choices.size() != sig2.choices.size()) {
        return true;  // different number of choices will always lead to a split
    }
    auto it2 = sig2.choices.begin();
    for (auto const& choice : sig1.choices) {
        if constexpr (Mode == SignatureMode::Exact) {
            if (choice != *it2) {
                return true;
            }
        } else {
            if (!choice.approxEqual(*it2, tolerance)) {
                return true;
            }
        }
        ++it2;
    }
    // Reaching this point means that we consider the signatures to be equal.
    return false;
}

template<typename ValueType, SignatureMode Mode>
void Signatures<ValueType, Mode>::extendQuotientData(QuotientData<ValueType>& quotientData) const {
    auto& signatureData = quotientData.signatureData.emplace();
    signatureData.toQuotientChoice.resize(model.getNumberOfChoices(), 0);
    signatureData.quotientChoiceGroupIndices.reserve(quotientData.toRepresentativeState.size() + 1);

    for (uint64_t quotientState = 0; quotientState < quotientData.toRepresentativeState.size(); ++quotientState) {
        auto const representativeState = quotientData.toRepresentativeState[quotientState];
        uint64_t const firstQuotientChoiceIndex = signatureData.toRepresentativeChoice.size();
        signatureData.quotientChoiceGroupIndices.push_back(firstQuotientChoiceIndex);
        // Relies on performSignatureBasedRefinement's postcondition that all cached signatures (including singleton blocks) are up to date.
        auto const& representativeSignature = stateSignatureCache[representativeState];

        std::vector<uint64_t> choiceSignatureToQuotientChoiceIndex(representativeSignature.choices.size(), std::numeric_limits<uint64_t>::max());
        for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(representativeState)) {
            auto [it, found] = representativeSignature.find(getChoiceSignature(choiceIndex), halfTolerance);
            STORM_LOG_ASSERT(found, "Expected to find the signature of representative state");
            uint64_t const choiceSignatureIndex = std::distance(representativeSignature.choices.begin(), it);
            if (choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex] == std::numeric_limits<uint64_t>::max()) {
                // We see this representative choice for the first time. Establish some mappings
                uint64_t const quotientChoiceIndex = signatureData.toRepresentativeChoice.size();
                choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex] = quotientChoiceIndex;
                signatureData.toRepresentativeChoice.push_back(choiceIndex);
                signatureData.toQuotientChoice[choiceIndex] = quotientChoiceIndex;
                // Fill distribution from signature
                auto& distr = signatureData.quotientChoiceDistributions.emplace_back();
                for (auto const& [successorBlock, value] : it->distr) {
                    distr.emplace(quotientData.toQuotientState[successorBlock.front()], value);
                }
            } else {
                // We have already seen this representative choice before.
                uint64_t const quotientChoiceIndex = choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex];
                signatureData.toQuotientChoice[choiceIndex] = quotientChoiceIndex;
            }
        }
        // Handle choice mappings for other blocks
        for (uint64_t const state : partition.getBlockOfElement(representativeState)) {
            if (state == representativeState) {
                continue;
            }
            for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(state)) {
                // This raw choice can be up to halfTolerance away from state's own canonical choice (established while updating state's signature during
                // refinement), which in turn can be up to halfTolerance away from the representative's canonical choice (established by SplitCondition).
                // Comparing the raw choice directly against the representative's canonical choice therefore spans two independent halfTolerance steps, i.e.
                // exactly the tolerance originally passed to this Signatures instance.
                auto [it, found] = representativeSignature.find(getChoiceSignature(choiceIndex), halfTolerance + halfTolerance);
                STORM_LOG_ASSERT(found, "Expected to find the signature of non-representative state");
                uint64_t const choiceSignatureIndex = std::distance(representativeSignature.choices.begin(), it);
                STORM_LOG_ASSERT(choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex] != std::numeric_limits<uint64_t>::max(),
                                 "Expected to have already seen this representative choice");
                uint64_t const quotientChoiceIndex = choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex];
                signatureData.toQuotientChoice[choiceIndex] = quotientChoiceIndex;
            }
        }
    }
    signatureData.quotientChoiceGroupIndices.push_back(signatureData.toRepresentativeChoice.size());
    signatureData.quotientChoiceGroupIndices.shrink_to_fit();
    signatureData.toRepresentativeChoice.shrink_to_fit();
    signatureData.quotientChoiceDistributions.shrink_to_fit();
}

template struct ChoiceSignatureCache<double>;
template struct ChoiceSignatureCache<storm::RationalNumber>;
template struct ChoiceSignatureCache<storm::RationalFunction>;

template class Signatures<double, SignatureMode::Exact>;
template class Signatures<double, SignatureMode::Approximative>;
template class Signatures<storm::RationalNumber, SignatureMode::Exact>;
template class Signatures<storm::RationalNumber, SignatureMode::Approximative>;
template class Signatures<storm::RationalFunction, SignatureMode::Exact>;

}  // namespace storm::bisimulation

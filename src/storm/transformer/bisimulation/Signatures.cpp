#include "storm/transformer/bisimulation/Signatures.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/ConstantsComparator.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::bisimulation {

template<typename ValueType, SignatureMode Mode>
std::strong_ordering StateSignature<ValueType, Mode>::ChoiceSignature::operator<=>(ChoiceSignature const& other) const {
    // It must hold that c_1 < c_2 < c_3 and c_1 ≈ c_3 implies c_1 ≈ c_2 ≈ c_3, where ≈ is the used "equality" for the Mode.
    // Since in all cases ≈ only holds if choiceClass and distr support coincide, we check those first
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
    // At this point only the distribution values could potentially differ.
    // In approximative mode, it is important that these are checked last so that similar distributions will also be close in the sorted choiceSignatures.
    // ValueType (e.g. RationalFunction, Interval) does not necessarily support <=>, so we fall back to </!= here.
    it2 = other.distr.begin();
    for (auto const& entry : distr) {
        if (entry.second != it2->second) {
            return entry.second < it2->second ? std::strong_ordering::less : std::strong_ordering::greater;
        }
        ++it2;
    }
    return std::strong_ordering::equal;
}

template<typename ValueType, SignatureMode Mode>
bool StateSignature<ValueType, Mode>::ChoiceSignature::operator==(ChoiceSignature const& other) const {
    return (*this <=> other) == std::strong_ordering::equal;
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
      tolerance(storm::utility::zero<ValueType>()),
      stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
Signatures<ValueType, Mode>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                        storm::bisimulation::Partition const& partition, ValueType const& tolerance)
    requires(Mode == SignatureMode::Approximative)
    : model(model), partition(partition), choiceClasses(choiceClasses), tolerance(tolerance), stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
void Signatures<ValueType, Mode>::updateStateSignature(uint64_t const stateIndex) {
    // Retrieves the mapping of blocks to (exact) aggregated probability
    using BlockDistributionType = typename Partition::OrderedBlockMap<ValueType>;
    auto getBlockDistribution = [this](uint64_t const choiceIndex) {
        BlockDistributionType distr;
        for (auto const& entry : model.getTransitionMatrix().getRow(choiceIndex)) {
            if (!storm::utility::isZero(entry.getValue())) {
                auto const emplace_res = distr.emplace(partition.getBlockOfElement(entry.getColumn()), entry.getValue());
                if (!emplace_res.second) {
                    emplace_res.first->second += entry.getValue();
                }
            }
        }
        return distr;
    };

    auto& choices = stateSignatureCache[stateIndex].choices;
    choices.clear();
    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        typename SignatureType::ChoiceSignature choiceSignature{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0,
                                                                .distr = getBlockDistribution(choiceIndex)};
        if constexpr (Mode == SignatureMode::Exact) {
            choices.emplace(std::move(choiceSignature), choiceIndex);
        } else {
            auto [it, inserted] = choices.emplace(std::move(choiceSignature), choiceIndex);
            if (inserted) {
                // Check if this distribution is approximately equal to one of the neighboring choices.
                // If yes, remove it again.
                if ((it != choices.begin() && it->first.approxEqual(std::prev(it)->first, tolerance)) ||
                    (std::next(it) != choices.end() && it->first.approxEqual(std::next(it)->first, tolerance))) {
                    choices.erase(it);
                }
            }
        }
    }
}

template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getSplitOrder() const -> SplitOrder {
    return SplitOrder{.signatures = stateSignatureCache};
}

template<typename ValueType, SignatureMode Mode>
bool Signatures<ValueType, Mode>::SplitOrder::operator()(uint64_t const state1, uint64_t const state2) const {
    SignatureType const& sig1 = signatures[state1];
    SignatureType const& sig2 = signatures[state2];
    // Do not compare the originating choice index (i.e. only compare the set of keys of the respective choices maps)
    if (auto const cmp = sig1.choices.size() <=> sig2.choices.size(); cmp != 0) {
        return cmp < 0;
    }
    auto it2 = sig2.choices.begin();
    for (auto const& [choice, _] : sig1.choices) {
        if (auto const cmp = choice <=> it2->first; cmp != 0) {
            return cmp < 0;
        }
        ++it2;
    }
    return false;
}
template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getSplitCondition() const -> SplitCondition {
    return SplitCondition{.signatures = stateSignatureCache, .tolerance = tolerance};
}

template<typename ValueType, SignatureMode Mode>
bool Signatures<ValueType, Mode>::SplitCondition::operator()(uint64_t const state1, uint64_t const state2) const {
    auto const& sig1 = signatures[state1];
    auto const& sig2 = signatures[state2];
    if (sig1.choices.size() != sig2.choices.size()) {
        return true;  // different number of choices will always lead to a split
    }
    auto it2 = sig2.choices.begin();
    for (auto const& [choice, _] : sig1.choices) {
        if constexpr (Mode == SignatureMode::Exact) {
            if (choice != it2->first) {
                return true;
            }
        } else {
            if (!choice.approxEqual(it2->first, tolerance)) {
                return true;
            }
        }
        ++it2;
    }
    // Reaching this point means that we consider the signatures to be equal.
    return false;
}

template class Signatures<double, SignatureMode::Exact>;
template class Signatures<double, SignatureMode::Approximative>;
template class Signatures<storm::RationalNumber, SignatureMode::Exact>;
template class Signatures<storm::RationalNumber, SignatureMode::Approximative>;
template class Signatures<storm::RationalFunction, SignatureMode::Exact>;

}  // namespace storm::bisimulation

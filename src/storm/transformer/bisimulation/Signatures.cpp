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
      tolerance(storm::utility::zero<ValueType>()),
      stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
Signatures<ValueType, Mode>::Signatures(storm::models::sparse::Model<ValueType> const& model, std::optional<std::vector<uint64_t>> const& choiceClasses,
                                        storm::bisimulation::Partition const& partition, ValueType const& tolerance)
    requires(Mode == SignatureMode::Approximative)
    : model(model), partition(partition), choiceClasses(choiceClasses), tolerance(tolerance), stateSignatureCache(model.getNumberOfStates()) {}

template<typename ValueType, SignatureMode Mode>
bool Signatures<ValueType, Mode>::ChoiceOrder::operator()(ChoiceSignature const& choice1, ChoiceSignature const& choice2) const {
    return choice1.compare(choice2, tolerance) == std::strong_ordering::less;
}

template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::mergeEquivalentChoices(uint64_t const stateIndex) const -> std::map<ChoiceSignature, std::vector<uint64_t>, ChoiceOrder> {
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

    std::map<ChoiceSignature, std::vector<uint64_t>, ChoiceOrder> groups(ChoiceOrder{.tolerance = tolerance});
    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        ChoiceSignature choiceSignature{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0, .distr = getBlockDistribution(choiceIndex)};
        auto [it, inserted] = groups.try_emplace(std::move(choiceSignature));
        if (!inserted) {
            it->second.push_back(choiceIndex);
        } else if constexpr (Mode == SignatureMode::Exact) {
            it->second.push_back(choiceIndex);
        } else {
            // Check if this distribution is approximately equal to one of the neighboring groups.
            // If yes, merge this choice into that group instead of keeping a separate one.
            if (it != groups.begin() && std::prev(it)->first.approxEqual(it->first, tolerance)) {
                std::prev(it)->second.push_back(choiceIndex);
                groups.erase(it);
            } else if (std::next(it) != groups.end() && std::next(it)->first.approxEqual(it->first, tolerance)) {
                std::next(it)->second.push_back(choiceIndex);
                groups.erase(it);
            } else {
                it->second.push_back(choiceIndex);
            }
        }
    }
    return groups;
}

template<typename ValueType, SignatureMode Mode>
void Signatures<ValueType, Mode>::updateStateSignature(uint64_t const stateIndex) {
    auto const mergedChoices = mergeEquivalentChoices(stateIndex);
    auto& choices = stateSignatureCache[stateIndex].choices;
    choices.clear();
    choices.reserve(mergedChoices.size());
    for (auto const& [signature, choiceIndices] : mergedChoices) {
        choices.push_back(signature);  // todo: check where we need the second entry ?
    }
}

template<typename ValueType, SignatureMode Mode>
auto Signatures<ValueType, Mode>::getSplitOrder() const -> SplitOrder {
    return SplitOrder{.signatures = stateSignatureCache, .tolerance = tolerance};
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

        // Merge the equivalent choices of the representative state; each resulting group of equivalent choices becomes exactly one quotient choice.
        auto const mergedRepresentativeChoices = mergeEquivalentChoices(representativeState);
        for (auto const& [signature, choiceIndices] : mergedRepresentativeChoices) {
            // Establish mappings between quotient choice and represented choice(s)
            uint64_t const quotientChoiceIndex = signatureData.toRepresentativeChoice.size();
            signatureData.toRepresentativeChoice.push_back(choiceIndices.front());
            for (uint64_t const choiceIndex : choiceIndices) {
                signatureData.toQuotientChoice[choiceIndex] = quotientChoiceIndex;
            }
            // Fill distribution from signature
            auto& distr = signatureData.quotientChoiceDistributions.emplace_back();
            for (auto const& [successorBlock, value] : signature.distr) {
                distr.emplace(quotientData.toQuotientState[successorBlock.front()], value);
            }
        }

        // Handle choice mappings for other blocks
        for (auto const state : partition.getBlockOfElement(representativeState)) {
            if (state == representativeState) {
                continue;
            }
            // All other states of the block must (by the invariant maintained during refinement) have the same number of distinct choices, in the same
        // (sorted) order.
            auto const mergedChoicesOfState = mergeEquivalentChoices(state);
            STORM_LOG_ASSERT(mergedChoicesOfState.size() == mergedRepresentativeChoices.size(),
                             "Expected states within the same block to have the same number of distinct (merged) choices.");
            STORM_LOG_ASSERT(([this, &mergedChoicesOfState, &mergedRepresentativeChoices]() {
                                 auto representativeIt = mergedChoicesOfState.begin();
                                 for (auto const& [signature, _] : mergedRepresentativeChoices) {
                                     if (signature.compare(representativeIt->first, tolerance) != 0) {
                                         return false;
                                     }
                                     ++representativeIt;
                                 }
                                 return true;
                             }()),
                             "Expected states within the same block to have the same (merged) choice signatures.");
            // Align mergeChoicesOfState with mergedRepresentativeChoices to obtain the quotient choice of each of the original choices.
            uint64_t quotientChoiceIndex = firstQuotientChoiceIndex;
            for (auto const& [_, choiceIndices] : mergedChoicesOfState) {
                for (uint64_t const choiceIndex : choiceIndices) {
                    signatureData.toQuotientChoice[choiceIndex] = quotientChoiceIndex;
                }
                ++quotientChoiceIndex;
            }
        }
    }
    signatureData.quotientChoiceGroupIndices.push_back(signatureData.toRepresentativeChoice.size());
    signatureData.quotientChoiceGroupIndices.shrink_to_fit();
    signatureData.toRepresentativeChoice.shrink_to_fit();
    signatureData.quotientChoiceDistributions.shrink_to_fit();
}

template class Signatures<double, SignatureMode::Exact>;
template class Signatures<double, SignatureMode::Approximative>;
template class Signatures<storm::RationalNumber, SignatureMode::Exact>;
template class Signatures<storm::RationalNumber, SignatureMode::Approximative>;
template class Signatures<storm::RationalFunction, SignatureMode::Exact>;

}  // namespace storm::bisimulation

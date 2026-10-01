#include "storm/transformer/bisimulation/Signatures.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <set>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/models/sparse/Model.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/utility/NumberTraits.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"
#include "storm/utility/matching.h"
#include "storm/utility/vector.h"

namespace storm::bisimulation {

namespace {
/*!
 * Shrinkens the intervals in the provided distribution to the feasible (aka coherent) parts of the intervals, i.e.,
 * those where the bounds can actually be instantiated to a valid probability distribution.
 */
template<typename IntervalType>
void makeFeasibleIntervalDistribution(std::span<std::pair<Partition::Block, IntervalType>> const distribution)
    requires IsIntervalType<IntervalType>
{
    using BaseType = storm::IntervalBaseType<IntervalType>;
    BaseType const one = storm::utility::one<BaseType>();

    // Helper function to set interval := interval cap [l,u]
    auto const intersect = [&one](IntervalType& interval, BaseType const& l, BaseType const& u) {
        interval.setLower(std::max(interval.lower(), l));
        interval.setUpper(std::min(interval.upper(), u));
        if constexpr (!storm::NumberTraits<BaseType>::IsExact) {
            // Ensure that the intersection is never empty, even if (for numerical reasons) the intervals are disjoint or if u < l.
            if (interval.lower() > interval.upper()) {
                interval = IntervalType((interval.lower() + interval.upper()) / (one + one));
            }
        }
        STORM_LOG_ASSERT(interval.lower() <= interval.upper(), "Restricting an interval to the coherent values made it empty.");
    };

    // A value can be at most one minus the sum of the lower bounds of the other values and at least one minus the sum of their upper bounds.

    // First, handle the fast cases with 1 or 2 entries.
    switch (distribution.size()) {
        case 0:
            return;
        case 1:
            distribution[0].second = one;
            return;
        case 2: {
            auto& v0 = distribution[0].second;
            auto& v1 = distribution[1].second;
            // Note that restricting v0 can only widen the restriction applied to v1 afterwards, so a single pass suffices.
            intersect(v0, one - v1.upper(), one - v1.lower());
            intersect(v1, one - v0.upper(), one - v0.lower());
            return;
        }
        default:
            break;  // continue with the general case below.
    }
    // Now the general case:
    IntervalType const sum =
        std::accumulate(distribution.begin(), distribution.end(), storm::utility::zero<IntervalType>(),
                        [](IntervalType const& acc, std::pair<Partition::Block, IntervalType> const& entry) { return acc + entry.second; });
    for (auto& entry : distribution) {
        auto& v = entry.second;
        BaseType const otherSumLower = sum.lower() - v.lower();
        BaseType const otherSumUpper = sum.upper() - v.upper();
        intersect(v, one - otherSumUpper, one - otherSumLower);
    }
}

template<typename ValueType>
std::strong_ordering compareValue(ValueType const& value1, ValueType const& value2) {
    if constexpr (storm::IsIntervalType<ValueType>) {
        auto cmpLower = compareValue(value1.lower(), value2.lower());
        if (cmpLower != std::strong_ordering::equal) {
            return cmpLower;
        }
        return compareValue(value1.upper(), value2.upper());
    } else {
        if (value1 < value2) {
            return std::strong_ordering::less;
        } else if (value1 == value2) {
            return std::strong_ordering::equal;
        } else {
            return std::strong_ordering::greater;
        }
    }
}
}  // namespace

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
Signatures<ValueType, Mode, QuotientValueType>::ChoiceSignatureCache::ChoiceSignatureCache(uint64_t const numStates)
    : values(numStates, storm::utility::zero<ValueType>()) {}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::ChoiceSignatureCache::addValue(Partition::Block const& block, ValueType const& value) {
    auto& current = values[block.front()];
    if (storm::utility::isZero(current)) {
        support.push_back(block);
    }
    current += value;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::ChoiceSignatureCache::extract(BlockDistributionView const dest) -> BlockDistributionView {
    STORM_LOG_ASSERT(support.size() <= dest.size(), "Destination span is too small to hold all entries.");

    // Prepare result: a sub-span of the destination view, which will be filled with the support.size() many (block, value) pairs.
    BlockDistributionView result = dest.first(support.size());

    // Sort the support for easier comparison of two distributions (i.e. choice signatures)
    std::sort(support.begin(), support.end(), Partition::BlockCompare());

    // Now write the (block, value) pairs while clearing the cache.
    {
        auto writeIt = result.begin();
        for (auto const& block : support) {
            auto& value = values[block.front()];
            if constexpr (std::is_same_v<ValueType, QuotientValueType>) {
                *writeIt = std::make_pair(block, std::move(value));
            } else {
                // The quotient abstracts the values of the model into intervals, so a single value only yields a point interval here.
                *writeIt = std::make_pair(block, storm::utility::convertNumber<QuotientValueType>(value));
            }
            ++writeIt;
            value = storm::utility::zero<ValueType>();  // reset value for next use of the cache.
        }
    }

    if (result.size() != dest.size()) {
        // Mark the end of the written entries with an empty block.
        // This is used to easily get the choice signature later without rebuilding it. See getChoiceSignature.
        dest[result.size()] = std::pair<Partition::Block, QuotientValueType>{};
    }
    support.clear();  // clear cache for next use.

    if constexpr (storm::IsIntervalType<QuotientValueType>) {
        makeFeasibleIntervalDistribution(result);
    }
    return result;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
std::strong_ordering Signatures<ValueType, Mode, QuotientValueType>::ConcreteChoiceSignature::compareStructure(ConcreteChoiceSignature const& other) const {
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
    return std::strong_ordering::equal;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::ConcreteChoiceSignature::compare(ConcreteChoiceSignature const& other) const -> ComparisonResult {
    if (auto const cmp = compareStructure(other); cmp != 0) {
        return cmp;
    }
    STORM_LOG_ASSERT(distr.size() == other.distr.size(), "The distributions of two signatures with the same structure must have the same size.");
    // Only the distribution values can still differ. QuotientValueType (e.g. RationalFunction, Interval) may lack <=>, so we fall back to </!= here.
    if constexpr (Mode == SignatureMode::Exact) {
        // Strong order: full lexicographic comparison.
        auto it2 = other.distr.begin();
        for (auto const& entry : distr) {
            if (auto const cmp = compareValue(entry.second, it2->second); cmp != std::strong_ordering::equal) {
                return cmp;
            }
            ++it2;
        }
    } else {
        // Weak order: compare only the first entry
        if (distr.empty()) {
            return ComparisonResult::equivalent;
        }
        QuotientValueType const& value = distr.front().second;
        QuotientValueType const& otherValue = other.distr.front().second;
        if constexpr (storm::IsIntervalType<QuotientValueType>) {
            // only compare the lower bound of the first entry
            if (auto cmp = compareValue(value.lower(), otherValue.lower()); cmp != std::strong_ordering::equal) {
                return cmp;
            }
        } else {
            if (auto cmp = compareValue(value, otherValue); cmp != std::strong_ordering::equal) {
                return cmp;
            }
        }
    }
    return ComparisonResult::equivalent;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
bool Signatures<ValueType, Mode, QuotientValueType>::ConcreteChoiceSignature::approximatelyEqual(ConcreteChoiceSignature const& other,
                                                                                                 ToleranceType const& tolerance) const
    requires(Mode == SignatureMode::Approximative)
{
    if (compareStructure(other) != std::strong_ordering::equal) {
        return false;
    }
    auto const near = [&tolerance](QuotientValueType const& value1, QuotientValueType const& value2) {
        return storm::utility::abs<QuotientValueType>(value1 - value2) <= tolerance;
    };
    auto otherIt = other.distr.begin();
    for (auto const& [block, value] : distr) {
        if (!near(value, otherIt->second)) {
            return false;
        }
        ++otherIt;
    }
    return true;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::AbstractChoiceSignature::enhance(AbstractChoiceSignature const& other)
    requires(Mode == SignatureMode::IntervalAbstraction)
{
    STORM_LOG_ASSERT(ConcreteChoiceSignature::compareStructure(other) == std::strong_ordering::equal,
                     "Cannot enhance two choice signatures with different structure.");
    using BaseType = storm::IntervalBaseType<QuotientValueType>;
    BaseType const zero = storm::utility::zero<BaseType>();

    // We keep track of how much the lower and upper bounds of the current signature have been increased and how much the other signature has been increased.
    // The latter is necessary since enhance is symmetric in the sense that a.enhance(b) and b.enhance(a) yield the same result.

    BaseType otherLowerDelta{other.lowerDelta}, otherUpperDelta{other.upperDelta};
    auto otherIt = other.distr.begin();
    for (auto& [block, value] : ConcreteChoiceSignature::distr) {
        // Enhance lower value
        BaseType const lowerDiff = otherIt->second.lower() - value.lower();
        if (lowerDiff < zero) {
            lowerDelta -= lowerDiff;
            value.setLower(otherIt->second.lower());
        } else {
            otherLowerDelta += lowerDiff;
        }
        // Enhance upper value
        BaseType const upperDiff = otherIt->second.upper() - value.upper();
        if (upperDiff > zero) {
            upperDelta += upperDiff;
            value.setUpper(otherIt->second.upper());
        } else {
            otherUpperDelta -= upperDiff;
        }
        ++otherIt;
    }
    lowerDelta = std::max<BaseType>(lowerDelta, otherLowerDelta);
    upperDelta = std::max<BaseType>(upperDelta, otherUpperDelta);
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
bool Signatures<ValueType, Mode, QuotientValueType>::AbstractChoiceSignature::isCompatibleWith(AbstractChoiceSignature const& other,
                                                                                               ToleranceType const& tolerance, bool requireContainsOther) const
    requires(Mode == SignatureMode::IntervalAbstraction)
{
    STORM_LOG_ASSERT(ConcreteChoiceSignature::compareStructure(other) == std::strong_ordering::equal, "Assuming same structure for compatible check.");

    using BaseType = storm::IntervalBaseType<QuotientValueType>;
    BaseType const zero = storm::utility::zero<BaseType>();

    BaseType thisLowerDelta{lowerDelta}, thisUpperDelta{upperDelta};
    BaseType otherLowerDelta{other.lowerDelta}, otherUpperDelta{other.upperDelta};
    auto otherIt = other.distr.begin();
    for (auto& [block, value] : ConcreteChoiceSignature::distr) {
        BaseType const lowerDiff = otherIt->second.lower() - value.lower();
        if (lowerDiff < zero) {
            if (requireContainsOther) {
                return false;
            }
            thisLowerDelta -= lowerDiff;
        } else {
            otherLowerDelta += lowerDiff;
        }
        BaseType const upperDiff = otherIt->second.upper() - value.upper();
        if (upperDiff > zero) {
            if (requireContainsOther) {
                return false;
            }
            thisUpperDelta += upperDiff;
        } else {
            otherUpperDelta -= upperDiff;
        }
        // Check if any of the deltas exceeds the tolerance
        if (thisLowerDelta > tolerance || thisUpperDelta > tolerance || otherLowerDelta > tolerance || otherUpperDelta > tolerance) {
            return false;
        }
        ++otherIt;
    }
    return true;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::StateSignature::find(ChoiceSignature const& signature, ToleranceType const& tolerance,
                                                                          [[maybe_unused]] bool requireContainsSignature) const
    -> std::pair<ChoiceSignatureIterator, bool> {
    if (choices.empty()) {
        return std::make_pair(choices.end(), false);
    }
    auto it = std::lower_bound(choices.begin(), choices.end(), signature, [](ChoiceSignature const& choice1, ChoiceSignature const& choice2) {
        return choice1.compare(choice2) == ChoiceSignature::ComparisonResult::less;
    });
    if constexpr (Mode == SignatureMode::Exact) {
        return std::make_pair(it, it != choices.end() && it->compare(signature) == std::strong_ordering::equal);
    } else {
        if (signature.distr.empty()) {
            // For empty choice signatures, it suffices to compare the structure.
            return std::make_pair(it, it != choices.end() && it->compareStructure(signature) == std::strong_ordering::equal);
        }
        return findWithHint(it, signature, tolerance, requireContainsSignature);
    }
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::StateSignature::findWithHint(ChoiceSignatureIterator hint, ChoiceSignature const& signature,
                                                                                  ToleranceType const& tolerance,
                                                                                  [[maybe_unused]] bool requireContainsSignature) const
    -> std::pair<ChoiceSignatureIterator, bool>
    requires(Mode != SignatureMode::Exact)
{
    // The input signature defines a (potentially empty) window inside the choices vector. That window is given by the entries that are compareStructure-equal
    // to signature and whose first entry is within tolerance of the first entry of signature.
    // As the choices are sorted by compare(), they can be devided into three consecutive, potentially empty segments:
    // the entries before the window, the window itself, and the entries after the window.

    STORM_LOG_ASSERT(hint >= choices.begin() && hint <= choices.end(), "Hint iterator must be within the range of choices.");
    STORM_LOG_ASSERT(!signature.distr.empty(), "Distributions must not be empty.");  // Don't deal with this special case here.

    using BaseType = storm::IntervalBaseType<QuotientValueType>;
    auto const near = [&tolerance](BaseType const& value1, BaseType const& value2) { return storm::utility::abs<BaseType>(value1 - value2) <= tolerance; };
    enum class Location { BeforeWindow, InWindow, AfterWindow };

    // Locate in which segment of the choices vector the given choiceIt is.
    auto const locate = [this, &signature, &near](ChoiceSignatureIterator choiceIt) {
        if (choiceIt == choices.end()) {
            return Location::AfterWindow;
        }
        auto const& choice = *choiceIt;
        if (auto const cmp = choice.compareStructure(signature); cmp != 0) {
            return cmp < 0 ? Location::BeforeWindow : Location::AfterWindow;
        }

        // Now the position only depends on the first entry as this is our sorting criterion for choice signatures (cf. ChoiceSignature::compare)
        // In case of intervals, we additionally just look at the lower bound of the first entry.
        auto getPosition = [&near](BaseType const& v, BaseType const& sigV) {
            if (near(v, sigV)) {
                return Location::InWindow;
            }
            return v < sigV ? Location::BeforeWindow : Location::AfterWindow;
        };

        auto const& value = choice.distr.front().second;
        auto const& signatureValue = signature.distr.front().second;
        if constexpr (Mode == SignatureMode::IntervalAbstraction) {
            return getPosition(value.lower(), signatureValue.lower());
        } else {
            return getPosition(value, signatureValue);
        }
    };

    // Returns true iff the choice is considered equivalent to the signature. Assumes that the given choice is in the window, in particular has the same
    // structure.
    auto const found = [&](ChoiceSignature const& choice) {
        if constexpr (Mode == SignatureMode::IntervalAbstraction) {
            return choice.isCompatibleWith(signature, tolerance, requireContainsSignature);
        } else {
            // The first entry is already assumed to be near, so we check the others.
            auto it1 = choice.distr.begin() + 1;
            auto it2 = signature.distr.begin() + 1;
            for (; it1 != choice.distr.end() && it2 != signature.distr.end(); ++it1, ++it2) {
                if (!near(it1->second, it2->second)) {
                    return false;
                }
            }
            return true;
        }
    };

    // Assumes an iterator that points inside the window and scans all window contents to the right (it excluded)
    auto const scanToRightBoundary = [&](ChoiceSignatureIterator it) {
        ++it;
        while (locate(it) == Location::InWindow) {
            if (found(*it)) {
                return std::make_pair(it, true);
            }
            ++it;
        }
        return std::make_pair(it, false);
    };
    // Assumes an iterator that points inside the window and scans all window contents to the left (it excluded)
    auto const scanToLeftBoundary = [&](ChoiceSignatureIterator it) {
        while (it != choices.begin()) {
            --it;
            if (locate(it) != Location::InWindow) {
                break;
            }
            if (found(*it)) {
                return std::make_pair(it, true);
            }
        }
        return std::make_pair(it, false);
    };

    // Callers usually pass a hint that lies in the window or is close to it. There are three cases for the hint:
    // a) inside, b) before, or c) after the window.
    auto const hintLocation = locate(hint);
    if (hintLocation == Location::InWindow) {  // case a), so we might need to scan in both directions
        if (found(*hint)) {
            return std::make_pair(hint, true);
        }
        if (auto resRight = scanToRightBoundary(hint); resRight.second) {
            return resRight;
        }
        if (auto resLeft = scanToLeftBoundary(hint); resLeft.second) {
            return resLeft;
        }
    } else {
        auto start = hint;
        auto startLocation = hintLocation;
        if (hintLocation == Location::BeforeWindow) {  // case b), find right boundary of window of the window
            while (startLocation == Location::BeforeWindow) {
                ++start;
                startLocation = locate(start);
            }
        } else {  // case c), find left boundary of the window
            while (startLocation == Location::AfterWindow && start != choices.begin()) {
                --start;
                startLocation = locate(start);
            }
        }
        if (startLocation == Location::InWindow) {
            if (found(*start)) {
                return std::make_pair(start, true);
            }
            if (auto result = hintLocation == Location::BeforeWindow ? scanToRightBoundary(start) : scanToLeftBoundary(start); result.second) {
                return result;
            }
        }
    }
    return std::make_pair(hint, false);
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::StateSignature::insert(ChoiceSignature const& choiceSignature, ToleranceType const& tolerance) {
    auto [it, found] = find(choiceSignature, tolerance);
    if (!found) {
        // Add fresh choice signature
        choices.insert(it, choiceSignature);
    } else if constexpr (Mode == SignatureMode::IntervalAbstraction) {
        // Suitable choice signature already exists. Enhance it!
        auto mutableIt = choices.begin() + std::distance(choices.cbegin(), it);
        mutableIt->enhance(choiceSignature);
    }
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
Signatures<ValueType, Mode, QuotientValueType>::Signatures(storm::models::sparse::Model<ValueType> const& model,
                                                           std::optional<std::vector<uint64_t>> const& choiceClasses,
                                                           storm::bisimulation::Partition const& partition)
    requires(Mode == SignatureMode::Exact)
    : choiceSignatureCache(partition.getNumberOfElements()),
      model(model),
      partition(partition),
      choiceClasses(choiceClasses),
      choiceDistributionStorage(model.getTransitionMatrix().getEntryCount()),
      halfTolerance(storm::utility::zero<ToleranceType>()),
      stateSignatureCache(model.getNumberOfStates()),
      tmpStateSignature(0) {}

namespace {
template<typename ValueType>
uint64_t getLargestRowGroupEntryCount(storm::storage::SparseMatrix<ValueType> const& matrix) {
    uint64_t maxSize = 0;
    for (uint64_t i = 0; i < matrix.getRowGroupCount(); ++i) {
        maxSize = std::max(maxSize, matrix.getRowGroup(i).getNumberOfEntries());
    }
    return maxSize;
}
}  // namespace

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
Signatures<ValueType, Mode, QuotientValueType>::Signatures(storm::models::sparse::Model<ValueType> const& model,
                                                           std::optional<std::vector<uint64_t>> const& choiceClasses,
                                                           storm::bisimulation::Partition const& partition, ToleranceType const& tolerance)
    requires(Mode != SignatureMode::Exact)
    : choiceSignatureCache(partition.getNumberOfElements()),
      model(model),
      partition(partition),
      choiceClasses(choiceClasses),
      choiceDistributionStorage(model.getTransitionMatrix().getEntryCount()),
      halfTolerance(tolerance / storm::utility::convertNumber<ToleranceType, uint64_t>(2)),
      stateSignatureCache(model.getNumberOfStates()),
      tmpStateSignature(Mode == SignatureMode::IntervalAbstraction ? getLargestRowGroupEntryCount(model.getTransitionMatrix()) : 0) {}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::buildChoiceSignature(uint64_t const choiceIndex) -> ChoiceSignature {
    auto const& matrix = model.getTransitionMatrix();
    // Use the choiceSignatureCache to compute the distribution over successor blocks
    for (auto const& entry : matrix.getRow(choiceIndex)) {
        if (!storm::utility::isZero(entry.getValue())) {
            choiceSignatureCache.addValue(partition.getBlockOfElement(entry.getColumn()), entry.getValue());
        }
    }
    // Identify the choice distribution storage.
    auto const row = matrix.getRow(choiceIndex);
    uint64_t const offset = std::distance(matrix.begin(), row.begin());
    BlockDistributionView distrStorage(choiceDistributionStorage.data() + offset, row.getNumberOfEntries());

    // Extract the resulting distribution from the cache.
    return ChoiceSignature{
        ConcreteChoiceSignature{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0, .distr = choiceSignatureCache.extract(distrStorage)}};
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::getChoiceSignature(uint64_t const choiceIndex) const -> ChoiceSignature {
    // Identify the choice distribution storage
    auto const& matrix = model.getTransitionMatrix();
    auto const row = matrix.getRow(choiceIndex);
    uint64_t const offset = std::distance(matrix.begin(), row.begin());
    BlockDistributionView distrStorage(choiceDistributionStorage.data() + offset, row.getNumberOfEntries());

    // Determine the right size by finding the first empty block (if any)
    uint64_t size = 0;
    for (auto const& [block, _] : distrStorage) {
        if (block.empty()) {
            break;
        }
        ++size;
    }

    return ChoiceSignature{ConcreteChoiceSignature{.choiceClass = choiceClasses ? (*choiceClasses)[choiceIndex] : 0, .distr = distrStorage.first(size)}};
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::updateStateSignature(uint64_t const stateIndex) {
    auto& sig = stateSignatureCache[stateIndex];
    sig.choices.clear();
    for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(stateIndex)) {
        sig.insert(buildChoiceSignature(choiceIndex), halfTolerance);
    }
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::getEquivalenceSplitOrder() const -> SplitOrder
    requires(Mode == SignatureMode::Exact)
{
    return SplitOrder(stateSignatureCache);
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::getStructuralSplitOrder() const -> SplitOrder
    requires(Mode != SignatureMode::Exact)
{
    return SplitOrder(stateSignatureCache);
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
bool Signatures<ValueType, Mode, QuotientValueType>::SplitOrder::operator()(uint64_t const state1, uint64_t const state2) const {
    auto const& sig1 = signatures[state1];
    auto const& sig2 = signatures[state2];
    if (auto const cmp = sig1.choices.size() <=> sig2.choices.size(); cmp != 0) {
        return cmp < 0;
    }
    auto it2 = sig2.choices.begin();
    for (auto const& choice : sig1.choices) {
        if constexpr (Mode == SignatureMode::Exact) {
            // In exact mode, do a full comparison of the choice signatures
            if (auto const cmp = choice.compare(*it2); cmp != 0) {
                return cmp < 0;
            }
        } else {
            // In approximative mode, only sort based on the structure of the choice signatures
            if (auto const cmp = choice.compareStructure(*it2); cmp != 0) {
                return cmp < 0;
            }
        }
        ++it2;
    }
    return false;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
auto Signatures<ValueType, Mode, QuotientValueType>::getClusteringSplitCondition() -> SplitCondition
    requires(Mode != SignatureMode::Exact)
{
    return SplitCondition(stateSignatureCache, halfTolerance, tmpStateSignature);
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::copyStructuralEquivalentStateSignature(uint64_t const& srcState, uint64_t const& dstState)
    requires(Mode == SignatureMode::IntervalAbstraction)
{
    auto const& src = stateSignatureCache[srcState];
    auto& dst = stateSignatureCache[dstState];
    STORM_LOG_ASSERT(src.choices.size() == dst.choices.size(), "Expected that source and destination have the same choice count.");
    for (uint64_t choiceIndex = 0; choiceIndex < src.choices.size(); ++choiceIndex) {
        auto const& srcChoice = src.choices[choiceIndex];
        auto& dstChoice = dst.choices[choiceIndex];
        STORM_LOG_ASSERT(srcChoice.compareStructure(dstChoice) == std::strong_ordering::equal,
                         "Expected that source and destination have the same choice structure.");
        std::copy(srcChoice.distr.begin(), srcChoice.distr.end(), dstChoice.distr.begin());
        dstChoice.lowerDelta = srcChoice.lowerDelta;
        dstChoice.upperDelta = srcChoice.upperDelta;
    }
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
bool Signatures<ValueType, Mode, QuotientValueType>::SplitCondition::operator()(uint64_t const anchorState, uint64_t const candidateState) const
    requires(Mode == SignatureMode::Approximative)
{
    auto const& sig1 = signatures[anchorState];
    auto const& sig2 = signatures[candidateState];
    // We can assume that sig1.choices and sig2.choices have the same pointwise structure (compareStructure-equivalent) as this comparison is only called after
    // splitting with respect to getStructuralSplitOrder.
    STORM_LOG_ASSERT(sig1.choices.size() == sig2.choices.size(), "SplitCondition should only be called for signatures with the same number of choices.");
    auto it1 = sig1.choices.begin();
    auto it2 = sig2.choices.begin();
    for (; it1 != sig1.choices.end() && it2 != sig2.choices.end(); ++it1, ++it2) {
        STORM_LOG_ASSERT(it1->distr.size() == it2->distr.size(), "SplitCondition should only be called for signatures with pointwise same choice structure.");
        if (it1->distr.empty()) {
            continue;
        }
        // Find a choice matching choice2 in sig1. As both signatures are sorted the same way, it1 is usually close to the window of choice2.
        auto const [choiceInSig1It, foundIn1] = sig1.findWithHint(it1, *it2, tolerance);
        if (!foundIn1) {
            return true;
        }
        // If it1 was not the matching choice, we have to check if it1 also has a matching choice in sig2
        if (it1 != choiceInSig1It && !sig2.findWithHint(it2, *it1, tolerance).second) {
            return true;
        }
    }
    // Reaching this point means that we consider the signatures to be equal, so no split.
    return false;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
bool Signatures<ValueType, Mode, QuotientValueType>::SplitCondition::operator()(uint64_t const anchorState, uint64_t const candidateState)
    requires(Mode == SignatureMode::IntervalAbstraction)
{
    // store the potential enhanced choice signature as a temporary, as we might not want to enhance the anchorState's signature if the candidateState does not
    // match.
    auto& sig1 = tmpStateSignature.load(signatures[anchorState]);
    auto const& sig2 = signatures[candidateState];
    // We can assume that sig1.choices and sig2.choices have the same pointwise structure (compareStructure-equivalent) as this comparison is only called after
    // splitting with respect to getStructuralSplitOrder.
    STORM_LOG_ASSERT(sig1.choices.size() == sig2.choices.size(), "SplitCondition should only be called for signatures with the same number of choices.");
    auto it1 = sig1.choices.begin();
    auto it2 = sig2.choices.begin();

    for (; it1 != sig1.choices.end() && it2 != sig2.choices.end(); ++it1, ++it2) {
        STORM_LOG_ASSERT(it1->distr.size() == it2->distr.size(), "SplitCondition should only be called for signatures with pointwise same choice structure.");
        if (it1->distr.empty()) {
            continue;
        }
        // Find a choice matching choice2 in sig1. As both signatures are sorted the same way, it1 is usually close to the window of choice2.
        auto const [choiceInSig1It, foundIn1] = sig1.findWithHint(it1, *it2, tolerance);
        if (foundIn1) {
            // findWithHint yields a const iterator, so we have to convert it into a mutable one to widen the found choice signature.
            auto mutableIt1 = sig1.choices.begin() + std::distance(sig1.choices.cbegin(), choiceInSig1It);
            mutableIt1->enhance(*it2);
        } else {
            // candidateState does not match with the anchorState: no match for choice2=*it2 in sig1
            return true;
        }
        // If it1 was not the matching choice, we have to check if it1 also has a matching choice in sig2
        if (it1 != choiceInSig1It) {
            auto const [choiceInSig2It, foundIn2] = sig2.findWithHint(it2, *it1, tolerance);
            if (foundIn2) {
                it1->enhance(*choiceInSig2It);
            } else {
                // candidateState does not match with the anchorState: no match for choice1=*it1 in sig2
                return true;
            }
        }
    }

    // Reaching this point means that we consider the signatures to be in the same block. Copy over the enhanced signature to the anchorState's signature.
    tmpStateSignature.store(signatures[anchorState]);
    return false;
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::addQuotientChoiceMapping(uint64_t const state, StateSignature const& representativeSignature,
                                                                              std::vector<uint64_t> const& choiceSignatureToQuotientChoiceIndex,
                                                                              std::vector<uint64_t>& toQuotientChoice) const {
    constexpr uint64_t Invalid = std::numeric_limits<uint64_t>::max();
    auto const& stateSignature = stateSignatureCache[state];
    STORM_LOG_ASSERT(stateSignature.choices.size() == representativeSignature.choices.size(), "States in a block have different number of choice signatures.");

    // Every choice of the state is represented by one of the state's own choice signatures. That representation is surjective, since each of those signatures
    // was established from a choice of this state. It thus suffices to find a bijective matching between the choice signatures of the state and those of the
    // representative: the composition of the two is surjective on the choice signatures of the representative as well. Such a surjective mapping is needed,
    // as otherwise a scheduler for the quotient model could not be translated back to a choice for this state.
    // In approximate mode, both steps are within halfTolerance, i.e. the values of a choice and those of the quotient choice representing it differ by at most
    // the tolerance this instance was created with.

    // Step 1: Find a matching between the choice signatures of the state and those of the representative.
    std::vector<uint64_t> choiceSignatureMatching;
    if constexpr (Mode == SignatureMode::Exact) {
        // The signatures of two states of the same block are equal entry-wise, so we can match them in order.
        choiceSignatureMatching = storm::utility::vector::buildVectorForRange<uint64_t>(0, stateSignature.choices.size());
        // Assert pointwise equality
        STORM_LOG_ASSERT(std::all_of(choiceSignatureMatching.begin(), choiceSignatureMatching.end(),
                                     [&stateSignature, &representativeSignature](uint64_t const i) {
                                         return stateSignature.choices[i].compare(representativeSignature.choices[i]) == std::strong_ordering::equal;
                                     }),
                         "In exact mode, choice signatures are expected to be equal for states in the same block.");
    } else if constexpr (Mode == SignatureMode::Approximative) {
        // In approximate mode, this is in fact a matching problem on a bipartite graph.
        auto optionalChoiceSignatureMatching = storm::utility::findPerfectMatching(
            stateSignature.choices.size(), [&stateSignature, &representativeSignature, this](uint64_t const stateChoice, uint64_t const representativeChoice) {
                return stateSignature.choices[stateChoice].approximatelyEqual(representativeSignature.choices[representativeChoice], halfTolerance);
            });
        // Such a matching does not have to exist: approximate equality is not transitive, so two choice signatures of the state can have the same, single
        // approximately equal partner in the representative's signature, which is not enough of a reason to split the block during refinement.
        STORM_LOG_THROW(optionalChoiceSignatureMatching.has_value(), storm::exceptions::UnexpectedException,
                        "Unable to match the choices of state " << state << " with the choices of the state that represents it in the quotient. "
                                                                << "Try again with a smaller tolerance.");
        choiceSignatureMatching = std::move(*optionalChoiceSignatureMatching);
    } else {
        // In interval abstraction mode, the signatures of two states of the same block are equal by construction as we have copied them over after enhancement.
        // Therefore, this is the same as in exact mode: we can match them in order.
        choiceSignatureMatching = storm::utility::vector::buildVectorForRange<uint64_t>(0, stateSignature.choices.size());
        // Assert pointwise equality by chacking compatibility with tolerance=0, which comes down to equality.
        STORM_LOG_ASSERT(std::all_of(choiceSignatureMatching.begin(), choiceSignatureMatching.end(),
                                     [&stateSignature, &representativeSignature](uint64_t const i) {
                                         auto const& c1 = stateSignature.choices[i];
                                         auto const& c2 = representativeSignature.choices[i];
                                         return (c1.compareStructure(c2) == std::strong_ordering::equal) &&
                                                c1.isCompatibleWith(c2, storm::utility::zero<ToleranceType>());
                                     }),
                         "In exact mode, choice signatures are expected to be equal for states in the same block.");
    }

    // Step 2: Map actual choice signatures to the representatives.
    auto const choiceIndices = model.getTransitionMatrix().getRowGroupIndices(state);
    for (uint64_t const choiceIndex : choiceIndices) {
        // Find the stored choice signature that represents the choice signature given by choiceIndex
        auto [it, found] = stateSignature.find(getChoiceSignature(choiceIndex), halfTolerance, true);
        STORM_LOG_ASSERT(found, "Expected to find the signature of a choice in the signature of its own state");
        uint64_t const indexInStateSignature = std::distance(stateSignature.choices.begin(), it);
        uint64_t const indexInRepresentativeSignature = choiceSignatureMatching[indexInStateSignature];
        STORM_LOG_ASSERT(choiceSignatureToQuotientChoiceIndex[indexInRepresentativeSignature] != Invalid,
                         "Expected to have already seen this representative choice");
        toQuotientChoice[choiceIndex] = choiceSignatureToQuotientChoiceIndex[indexInRepresentativeSignature];
    }

    // Finally, assert that we indeed have a surjective mapping
    STORM_LOG_ASSERT(([&]() {
                         std::set<uint64_t> seenQuotientChoices;
                         for (uint64_t const choiceIndex : choiceIndices) {
                             seenQuotientChoices.insert(toQuotientChoice[choiceIndex]);
                         }
                         return seenQuotientChoices.size() == representativeSignature.choices.size();
                     }()),
                     "The mapping from state choices to quotient choices is not surjective.");
}

template<typename ValueType, SignatureMode Mode, typename QuotientValueType>
void Signatures<ValueType, Mode, QuotientValueType>::extendQuotientData(QuotientData<QuotientValueType>& quotientData,
                                                                        bool const createQuotientChoiceMapping) const {
    auto& signatureData = quotientData.signatureData.emplace();
    if (createQuotientChoiceMapping) {
        quotientData.toQuotientChoice.emplace(model.getNumberOfChoices(), 0);
    }
    signatureData.quotientChoiceGroupIndices.reserve(quotientData.toRepresentativeState.size() + 1);

    // Scratch space for mappings, reused for every state.
    std::vector<uint64_t> choiceSignatureToQuotientChoiceIndex;  // Maps choiceSignature indices to quotient model choice indices
    constexpr uint64_t Invalid = std::numeric_limits<uint64_t>::max();

    for (uint64_t quotientState = 0; quotientState < quotientData.toRepresentativeState.size(); ++quotientState) {
        auto const representativeState = quotientData.toRepresentativeState[quotientState];
        uint64_t const firstQuotientChoiceIndex = signatureData.toRepresentativeChoice.size();
        signatureData.quotientChoiceGroupIndices.push_back(firstQuotientChoiceIndex);
        // Relies on performSignatureBasedRefinement's postcondition that all cached signatures (including singleton blocks) are up to date.
        auto const& representativeSignature = stateSignatureCache[representativeState];

        // At the quotientState, the choices should appear in the same order as the represented choices at the representativeState. This prevents unexpected
        // effects downstream (e.g. different quality of the initial policy in policy iteration). For Markov automata with hybrid states, this is even necessary
        // to ensure that the Markovian choice remains the first one.
        // The choiceSignatureToQuotientChoiceIndex mapping below permutes the choices as they appear in representativeSignature into the right order.
        choiceSignatureToQuotientChoiceIndex.assign(representativeSignature.choices.size(), Invalid);
        for (uint64_t const choiceIndex : model.getTransitionMatrix().getRowGroupIndices(representativeState)) {
            // Find the stored choice signature that represents the choice signature given by choiceIndex
            auto [it, found] = representativeSignature.find(getChoiceSignature(choiceIndex), halfTolerance, true);
            STORM_LOG_ASSERT(found, "Expected to find the signature of representative state");
            uint64_t const choiceSignatureIndex = std::distance(representativeSignature.choices.begin(), it);
            if (choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex] == Invalid) {
                // We see this representative choice for the first time. Establish some mappings
                uint64_t const quotientChoiceIndex = signatureData.toRepresentativeChoice.size();
                choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex] = quotientChoiceIndex;
                signatureData.toRepresentativeChoice.push_back(choiceIndex);
                if (createQuotientChoiceMapping) {
                    (*quotientData.toQuotientChoice)[choiceIndex] = quotientChoiceIndex;
                }
                auto& distr = signatureData.quotientChoiceDistributions.emplace_back();
                for (auto const& [successorBlock, value] : it->distr) {
                    distr.emplace(quotientData.toQuotientState[successorBlock.front()], value);
                }
            } else {
                // We have already seen this representative choice before.
                if (createQuotientChoiceMapping) {
                    (*quotientData.toQuotientChoice)[choiceIndex] = choiceSignatureToQuotientChoiceIndex[choiceSignatureIndex];
                }
            }
        }
        if (createQuotientChoiceMapping) {
            // Handle choice mappings for the other states of the block. This is comparatively expensive, since (unlike the representative's own choices
            // above) it visits every choice of every non-representative state, so it is skipped entirely if the caller does not need toQuotientChoice.
            for (uint64_t const state : partition.getBlockOfElement(representativeState)) {
                if (state != representativeState) {
                    addQuotientChoiceMapping(state, representativeSignature, choiceSignatureToQuotientChoiceIndex, *quotientData.toQuotientChoice);
                }
            }
        }
    }
    signatureData.quotientChoiceGroupIndices.push_back(signatureData.toRepresentativeChoice.size());
    signatureData.quotientChoiceGroupIndices.shrink_to_fit();
    signatureData.toRepresentativeChoice.shrink_to_fit();
    signatureData.quotientChoiceDistributions.shrink_to_fit();
}

// Explicit instantiations for QuotientValueType == ValueType
template class Signatures<double, SignatureMode::Exact>;
template class Signatures<double, SignatureMode::Approximative>;
template class Signatures<storm::RationalNumber, SignatureMode::Exact>;
template class Signatures<storm::RationalNumber, SignatureMode::Approximative>;
template class Signatures<storm::RationalFunction, SignatureMode::Exact>;
template class Signatures<storm::Interval, SignatureMode::Exact>;
template class Signatures<storm::Interval, SignatureMode::IntervalAbstraction>;
template class Signatures<storm::RationalInterval, SignatureMode::Exact>;
template class Signatures<storm::RationalInterval, SignatureMode::IntervalAbstraction>;

// Explicit instantiations for QuotientValueType == IntervalType<ValueType> (for interval abstraction)
// Exact mode is not meaningful in this case, as that would mean that we are never allowed to abstract values into intervals.
template class Signatures<double, SignatureMode::IntervalAbstraction, storm::Interval>;
template class Signatures<storm::RationalNumber, SignatureMode::IntervalAbstraction, storm::RationalInterval>;

}  // namespace storm::bisimulation

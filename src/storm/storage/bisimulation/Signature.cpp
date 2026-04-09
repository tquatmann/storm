#include "Signature.h"

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace storage {
namespace bisimulation {

template<typename ValueType>
void Signature<ValueType>::addBlockProbability(size_t blockId, ValueType probability) {
    auto it = blockProbabilities.find(blockId);
    if (it != blockProbabilities.end()) {
        it->second += probability;
    } else {
        blockProbabilities.emplace(blockId, probability);
    }
}

template<typename ValueType>
size_t Signature<ValueType>::computeHash() const {
    std::size_t seed = 0;
    // if (storm::IsIntervalType<ValueType>) {
    //     for (const auto& [blockId, prob] : blockProbabilities) {
    //         boost::hash_combine(seed, blockId);
    //         boost::hash_combine(seed, prob);
    //     }
    // }

    return seed;
}

template<typename ValueType>
std::string Signature<ValueType>::toString() const {
    std::string result;
    for (const auto& [blockId, prob] : blockProbabilities) {
        if constexpr (std::is_same_v<decltype(prob), const double>) {
            result += std::to_string(blockId) + ":" + std::to_string(prob) + ",";
        }
    }
    return result;
}

template class Signature<double>;
template class Signature<storm::Interval>;
template class Signature<storm::RationalNumber>;
template class Signature<storm::RationalFunction>;

}  // namespace bisimulation
}  // namespace storage
}  // namespace storm
#pragma once

#include <boost/container/flat_map.hpp>
#include <iostream>

#include "storm/adapters/IntervalAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace storage {
namespace bisimulation {

template<typename ValueType>
class Signature {
   public:
    boost::container::flat_map<std::size_t, ValueType> blockProbabilities;

    Signature() = default;

    bool operator==(const Signature &other) const {
        return blockProbabilities == other.blockProbabilities;
    }

    bool operator<(const Signature &other) const {
        return blockProbabilities < other.blockProbabilities;
    }

    void addBlockProbability(size_t blockId, ValueType probability);

    size_t computeHash() const;

    std::string toString() const;
};
}  // namespace bisimulation
}  // namespace storage
}  // namespace storm

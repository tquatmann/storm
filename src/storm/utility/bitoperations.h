#pragma once
#include <array>
#include <bit>
#include <cstdint>
#include <ranges>

#include "storm/utility/macros.h"

namespace storm::utility {

/*!
 * Swaps the byte representation of the given value, (i.e., swaps endianness)
 * Taken from https://en.cppreference.com/w/cpp/numeric/byteswap
 * Can be replaced once we are at C++23
 */
template<typename T>
inline T byteSwap(T const t) {
    // Note: c++23's std::bit_cast does not allow floating point types. Potentially because they don't have a unique object representation (NaN).
    static_assert(std::has_unique_object_representations_v<T> || std::is_same_v<T, double> || std::is_same_v<T, float>, "T may not have padding bits");
    auto value_representation = std::bit_cast<std::array<std::byte, sizeof(T)>>(t);
    std::ranges::reverse(value_representation);
    return std::bit_cast<T>(value_representation);
}

}  // namespace storm::utility

/**
 * \return 2^n - 1
 */
inline uint64_t smallestIntWithNBitsSet(uint64_t n) {
    STORM_LOG_ASSERT(n < 64, "Input is too large.");
    if (n == 0)
        return 0;
    return (1ul << n) - 1;
}

/**
 * The next bit permutation in a lexicographical sense.
 *
 * Example: 00010011, 00010101, 00010110, 00011001,
 * 00011010, 00011100, 00100011, and so forth
 *
 * From https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
 */
inline uint64_t nextBitPermutation(uint64_t v) {
    if (v == 0)
        return 0;
    uint64_t t = (v | (v - 1)) + 1;
    return t | ((((t & -t) / (v & -v)) >> 1) - 1);
}

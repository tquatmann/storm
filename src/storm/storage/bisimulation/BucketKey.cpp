#include "BucketKey.h"
#include <cmath>
#include <algorithm>
#include "storm/storage/bisimulation/Partition.h" // defines Partition, forEachBlock, etc.

namespace storm {
namespace storage {
namespace bisimulation {

bool BucketKey::operator<(BucketKey const& o) const noexcept {
    if (entries.size() != o.entries.size())
        return entries.size() < o.entries.size();
    for (std::size_t i = 0; i < entries.size(); ++i) {
        if (entries[i] < o.entries[i])
            return true;
        if (o.entries[i] < entries[i])
            return false;
    }
    return false;
}

bool BucketKey::operator==(BucketKey const& o) const noexcept {
    return entries == o.entries;
}

std::size_t BucketKey::Hash::operator()(BucketKey const& k) const noexcept {
    std::size_t h = 1469598103934665603ULL;  // FNV-ish
    auto mix = [&](std::size_t x) { h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); };
    for (auto const& e : k.entries) {
        mix(static_cast<std::size_t>(e.rep));
        mix(static_cast<std::size_t>(e.code.lo));
        mix(static_cast<std::size_t>(e.code.hi));
    }
    return h;
}

double clamp01(double x) noexcept {
    if (x < 0.0)
        return 0.0;
    if (x > 1.0)
        return 1.0;
    return x;
}

int binIndex(double x, double delta) noexcept {
    // defend against tiny FP noise on boundaries
    return static_cast<int>(std::floor(clamp01(x) / delta + 1e-12));
}

BucketCode makeBucket(double L, double U, double delta) noexcept {
    return BucketCode{binIndex(L, delta), binIndex(U, delta)};
}

BucketKeyBuilder::BucketKeyBuilder(Partition const& part, double delta, ProjectionFn proj) : partition_(part), delta_(delta), proj_(std::move(proj)) {}

BucketKey BucketKeyBuilder::buildForState(std::size_t s) const {
    BucketKey key;
    key.entries.reserve(partition_.getNumberOfBlocks());

    // Collect (rep, code) pairs in whatever order...
    partition_.forEachBlock([&](std::span<const std::uint64_t> block) {
        std::uint64_t rep = block.front();  // canonical ID of block
        auto [L, U] = proj_(s, block);
        key.entries.push_back({rep, makeBucket(L, U, delta_)});
    });

    // ...then sort by rep to canonicalize order.
    std::sort(key.entries.begin(), key.entries.end(), [](auto const& a, auto const& b) { return a.rep < b.rep || (a.rep == b.rep && a.code < b.code); });

    return key;
}

}  // namespace bisimulation
}  // namespace storage
}  // namespace storm
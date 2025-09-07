#pragma once
#include <vector>
#include <functional>
#include <cstddef>
#include <span>
#include <cstdint>

namespace storm {
namespace storage {
namespace bisimulation {

// A single bucket code for one block: (lower_idx, upper_idx)
struct BucketCode {
    int lo = 0;
    int hi = 0;

    bool operator<(BucketCode const& o) const noexcept {
        return (lo < o.lo) || (lo == o.lo && hi < o.hi);
    }
    bool operator==(BucketCode const& o) const noexcept {
        return lo == o.lo && hi == o.hi;
    }
};

// We store (rep-id, code) for each block and compare lexicographically by rep-id then code.
// This makes keys stable w.r.t. arbitrary block enumeration orders.
struct BucketKey {
    struct Entry {
        std::uint64_t rep = 0; // canonical representative of the block (e.g., block.front())
        BucketCode code{};
        bool operator<(Entry const& o) const noexcept {
            if (rep != o.rep) return rep < o.rep;
            return code < o.code;
        }
        bool operator==(Entry const& o) const noexcept {
            return rep == o.rep && code == o.code;
        }
    };

    std::vector<Entry> entries;

    bool operator<(BucketKey const& o) const noexcept;
    bool operator==(BucketKey const& o) const noexcept;

    struct Hash {
        std::size_t operator()(BucketKey const& k) const noexcept;
    };
};

// Helpers to build bucket codes
double clamp01(double x) noexcept;
int    binIndex(double x, double delta) noexcept;
BucketCode makeBucket(double L, double U, double delta) noexcept;

// Forward decl: we only need the type for the builder API.
class Partition;

// Builder that knows how to iterate the current partition and compute a key.
class BucketKeyBuilder {
   public:
    using ProjectionFn =
        std::function<
            // returns [L,U] as a pair<double,double> for state s to block C
            std::pair<double,double>(std::size_t /*state*/, std::span<const std::uint64_t> /*block*/)
            >;

    BucketKeyBuilder(Partition const& part, double delta, ProjectionFn proj);

    BucketKey buildForState(std::size_t s) const;

   private:
    Partition const& partition_;
    double delta_;
    ProjectionFn proj_;
};

}  // namespace bisimulation
}  // namespace storage
}  // namespace storm

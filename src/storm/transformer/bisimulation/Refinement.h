#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/transformer/bisimulation/Partition.h"

namespace storm::bisimulation {

// idea for signature-based refinement:
// flag blocks that are known to be stable w.r.t. the current partition
// Invariant: all unflaged blocks are in the queue (there might also be flagged ones in the queue)
// take the smallest block B from the queue
// Split B  into B=B_1 cup B_2 cup ... cup B_n using signature (if B is flagged, we already know that B=B_n=B_1)
// If refinement of B had no effect, i.e. it  was already stable (B_n=B_1=B) then flag B
// Compute for each predecessor state s of B an uint64_t s_i with the following interpretation:
// s_i == 0: s has no transition to B
// 0 < s_i <= n: s_i has a transition to B_i and no other sub-block of B
// s_i > n: s_i is a hash encoding the subset of {B_1,...,B_n} that s has transitions to (can enforce this by always setting the most significant bit to true)
// Note: if B is flagged, then s_i is either 0 or 1
// Then, for each predecessor block B' of B:
//   if B' and B are flagged: continue (assert that  s_i=1 for all s in B')
//   split B' according to s_i.
//   if B' is flagged (and B is not flagged):
//      - there shouldn't be any s in B' with s_i = 0 (assert that! could also make the split of B' above more efficient)
//      - for subblocks of B' with 0 < s_i <= n, all transitions that lead to B also lead to B_i so the split of B does not introduce any distinguishing
//           behavior within such blocks. flag those subblocks of B' as well
//      - subblocks of B' with s_i > n should no longer be flagged (including the case where B' was not split)
//    if B' was split above: add all subblocks of B' to the queue
//    else if B' is unflagged: if it was already unflagged, assert that it's already in the queue. Add B' to the queue
//   Nothing needs to be done if B' is flagged and not split. It might or might not be already in the queue

template<typename ValueType>
void performPartitionRefinement(storm::models::sparse::Model<ValueType> const& model, storm::bisimulation::Partition& partition);

}  // namespace storm::bisimulation

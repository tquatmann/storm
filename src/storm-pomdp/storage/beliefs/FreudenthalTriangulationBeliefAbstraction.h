#pragma once

#include <vector>

#include "storm-pomdp/storage/beliefs/BeliefBuilder.h"
#include "storm-pomdp/storage/beliefs/BeliefNumerics.h"
#include "storm-pomdp/storage/beliefs/BeliefTypes.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm::pomdp::beliefs {

template<typename BeliefType>
class FreudenthalTriangulationBeliefAbstraction {
   public:
    using BeliefValueType = typename BeliefType::ValueType;

    enum class Mode { Static, Dynamic };

    FreudenthalTriangulationBeliefAbstraction(std::vector<BeliefValueType> const& observationTriangulationResolutions, Mode mode) : mode(mode) {
        observationResolutions.reserve(observationResolutions.size());
        for (auto const& res : observationTriangulationResolutions) {
            observationResolutions.emplace_back(storm::utility::ceil<BeliefValueType>(res));
            STORM_LOG_ASSERT(observationResolutions.back() > 0, "Expected that the resolution is a positive integer. Got " << res << " instead.");
        }
    }

    template<typename AbstractCallback>
    void abstract(BeliefValueType&& probabilityFactor, BeliefType&& belief, AbstractCallback const& callback) {
        // Quickly triangulate Dirac beliefs
        if (belief.size() == 1u) {
            callback(std::move(probabilityFactor), std::move(belief));
        } else {
            switch (mode) {
                case Mode::Static:
                    abstractStatic(probabilityFactor, belief, callback, observationResolutions[belief.observation()]);
                    break;
                case Mode::Dynamic:
                    abstractDynamic(probabilityFactor, belief, callback);
                    break;
                default:
                    STORM_LOG_ASSERT(false, "Invalid triangulation mode.");
            }
        }
    }

   private:
    template<typename AbstractCallback>
    void abstractStatic(BeliefValueType const& probabilityFactor, BeliefType const& belief, AbstractCallback const& callback,
                        BeliefValueType const& staticResolution) {
        struct FreudenthalDiff {
            BeliefValueType diff;       // d[i]
            BeliefStateType dimension;  // i
            bool operator>(FreudenthalDiff const& other) const {
                if (diff != other.diff) {
                    return diff > other.diff;
                } else {
                    return dimension < other.dimension;
                }
            }
        };

        STORM_LOG_ASSERT(storm::utility::isInteger(staticResolution), "Expected an integer resolution");
        STORM_LOG_ASSERT(staticResolution > 0, "Expected a positive resolution");
        BeliefStateType numEntries = belief.size();
        // This is the Freudenthal Triangulation as described in Lovejoy (a whole lotta math)
        // Probabilities will be triangulated to values in 0/N, 1/N, 2/N, ..., N/N
        // Variable names are mostly based on the paper
        // However, we speed this up a little by exploiting that belief states usually have sparse support (i.e. numEntries is much smaller than
        // pomdp.getNumberOfStates()). Initialize diffs and the first row of the 'qs' matrix (aka v)
        std::set<FreudenthalDiff, std::greater<>> sorted_diffs;  // d (and p?) in the paper
        std::vector<BeliefValueType> qsRow;                      // Row of the 'qs' matrix from the paper (initially corresponds to v
        qsRow.reserve(numEntries);
        std::vector<BeliefStateType> toOriginalIndicesMap;  // Maps 'local' indices to the original pomdp state indices
        toOriginalIndicesMap.reserve(numEntries);
        BeliefValueType x = staticResolution;
        for (auto const& entry : belief) {
            qsRow.push_back(storm::utility::floor(x));                              // v
            sorted_diffs.emplace({x - qsRow.back(), toOriginalIndicesMap.size()});  // x-v
            toOriginalIndicesMap.push_back(entry.first);
            x -= entry.second * staticResolution;
        }
        // Insert a dummy 0 column in the qs matrix so the loops below are a bit simpler
        qsRow.push_back(storm::utility::zero<BeliefValueType>());

        auto currentSortedDiff = sorted_diffs.begin();
        auto previousSortedDiff = sorted_diffs.end();
        --previousSortedDiff;
        for (BeliefStateType i = 0; i < numEntries; ++i) {
            // Compute the weight for the grid points
            BeliefValueType weight = previousSortedDiff->diff - currentSortedDiff->diff;
            if (i == 0) {
                // The first weight is a bit different
                weight += storm::utility::one<BeliefValueType>();
            } else {
                // 'compute' the next row of the qs matrix
                qsRow[previousSortedDiff->dimension] += storm::utility::one<BeliefValueType>();
            }
            if (!BeliefNumerics<BeliefValueType>::isZero(weight)) {
                // build the grid point
                BeliefBuilder<BeliefType> builder;
                builder.reserve(numEntries);
                builder.setObservation(belief.observation());
                for (BeliefStateType j = 0; j < numEntries; ++j) {
                    BeliefValueType gridPointEntry = qsRow[j] - qsRow[j + 1];
                    if (!BeliefNumerics<BeliefValueType>::isZero(gridPointEntry)) {
                        builder.addValue(toOriginalIndicesMap[j], gridPointEntry / staticResolution);
                    }
                }
                callback(static_cast<BeliefValueType>(weight * probabilityFactor), builder.build());
            }
            previousSortedDiff = currentSortedDiff++;
        }
    }

    template<typename AbstractCallback>
    void abstractDynamic(BeliefValueType const& probabilityFactor, BeliefType const& belief, AbstractCallback const& callback) {
        // Find the best resolution for this belief, i.e., N such that the largest distance between one of the belief values to a value in {i/N | 0 ≤ i ≤ N} is
        // minimal
        auto const resolution = observationResolutions[belief.observation()];
        STORM_LOG_ASSERT(storm::utility::isInteger(resolution), "Expected an integer resolution");
        BeliefValueType finalResolution = resolution;
        BeliefValueType finalResolutionDist = storm::utility::one<BeliefValueType>();
        // We don't need to check resolutions that are smaller than the maximal resolution divided by 2 as we already checked multiples of these
        for (BeliefValueType currResolution = resolution; currResolution > resolution / 2; --currResolution) {
            BeliefValueType currDist = storm::utility::zero<BeliefValueType>();
            bool const newBest = belief.allOf([&currDist, &currResolution](BeliefStateType const&, BeliefValueType const& val) {
                currDist += storm::utility::abs<BeliefValueType>(val - storm::utility::round<BeliefValueType>(val * currResolution) / currResolution);
                return currDist <= finalResolutionDist;  // continue as long as the current dist is still smaller than the smallest dist
            });
            if (newBest) {
                STORM_LOG_ASSERT(currDist <= finalResolutionDist, "Expected a smaller resolution");
                finalResolution = currResolution;
                finalResolutionDist = currDist;
                if (BeliefNumerics<BeliefValueType>::isZero(finalResolutionDist)) {
                    break;
                }
            }
        }
        abstractStatic(probabilityFactor, belief, callback, finalResolution);
    }

   private:
    std::vector<BeliefValueType> observationResolutions;
    Mode const mode;
};
}  // namespace storm::pomdp::beliefs
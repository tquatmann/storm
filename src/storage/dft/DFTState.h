#ifndef DFTSTATE_H
#define	DFTSTATE_H

#include "../BitVector.h"
#include "DFTElementState.h"

#include <sstream>

namespace storm {
    namespace storage {
       
        class DFT; 
        template<typename T>
        class DFTBE;
        
        
        
        class DFTState {
            friend struct std::hash<DFTState>;
        private:
            storm::storage::BitVector mStatus;
            size_t mId;
            std::vector<size_t> mInactiveSpares;
            std::vector<size_t> mIsCurrentlyFailableBE;
            std::vector<size_t> mUsedRepresentants;
            bool mValid = true;
            const DFT& mDft;
            
        public:
            DFTState(DFT const& dft, size_t id);
            
            DFTElementState getElementState(size_t id) const;

            int getElementStateInt(size_t id) const;

            size_t getId() const;

            void setId(size_t id);
            
            bool isOperational(size_t id) const;
            
            bool hasFailed(size_t id) const;
            
            bool isFailsafe(size_t id) const ;
            
            bool dontCare(size_t id) const;
            
            void setFailed(size_t id);
            
            void setFailsafe(size_t id);
            
            void setDontCare(size_t id);
            
            void beNoLongerFailable(size_t id);
            
            void activate(size_t repr);
            
            bool isActiveSpare(size_t id) const;
            
            void markAsInvalid() {
                mValid = false;
            }
           
            bool isInvalid() {
                return !mValid;
            }
            
            storm::storage::BitVector const& status() const {
                return mStatus;
            }
           
            /**
             * This method gets the usage information for a spare
             * @param id Id of the spare
             * @return The child that currently is used.
             */
            uint_fast64_t uses(size_t id) const;
            
            /**
             * This method is commonly used to get the usage information for spares. 
             * @param from Starting index where the usage info is.
             * @return The child that currently is used.
             */
            uint_fast64_t extractUses(size_t from) const;
            
            /**
             * Checks whether an element is currently used.
             * @param child The id of the child for which we want to know whether it is currently used.
             * @return true iff it is currently used by any of the spares.
             */
            bool isUsed(size_t child);
            
            /**
             * Sets to to the usageIndex which child is now used.
             * @param usageIndex 
             * @param child
             */
            void setUsesAtPosition(size_t usageIndex, size_t child);
            
            
            
            bool claimNew(size_t spareId, size_t usageIndex, size_t currentlyUses, std::vector<size_t> const& childIds);
            
            bool hasOutgoingEdges() const {
                return !mIsCurrentlyFailableBE.empty();
            }
            
            size_t nrFailableBEs() const {
                return mIsCurrentlyFailableBE.size();
            }
            
            std::pair<std::shared_ptr<DFTBE<double>>, bool> letNextBEFail(size_t smallestIndex = 0);
            
            std::string getCurrentlyFailableString() {
                std::stringstream stream;
                auto it = mIsCurrentlyFailableBE.begin();
                stream << "{";
                if(it != mIsCurrentlyFailableBE.end()) {
                    stream << *it;
                }
                ++it;
                while(it != mIsCurrentlyFailableBE.end()) {
                    stream << ", " << *it;
                    ++it;
                }
                stream << "}";
                return stream.str();
            }

            friend bool operator==(DFTState const& a, DFTState const& b) {
                return a.mStatus == b.mStatus;
            }
        };

    }
}

namespace std {
    template<>
    struct hash<storm::storage::DFTState> {
        size_t operator()(storm::storage::DFTState const& s) const {
            return hash<storm::storage::BitVector>()(s.mStatus);
        }
    };
}

#endif	/* DFTSTATE_H */


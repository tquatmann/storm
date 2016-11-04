#pragma once

#include "src/builder/jit/Distribution.h"

namespace storm {
    namespace builder {
        namespace jit {
            
            template <typename IndexType, typename ValueType>
            class Choice {
            public:
                Choice();
                
                void setMarkovian(bool value);
                bool isMarkovian() const;
                
                void add(DistributionEntry<IndexType, ValueType> const& entry);
                void add(IndexType const& index, ValueType const& value);
                void add(Choice<IndexType, ValueType>&& choice);
                
                Distribution<IndexType, ValueType> const& getDistribution() const;
                void divideDistribution(ValueType const& value);
                
                /*!
                 * Adds the given value to the reward associated with this choice.
                 */
                void addReward(ValueType const& value);
                
                /*!
                 * Adds the given value to the reward with the given index of this choice.
                 */
                void addReward(uint64_t index, ValueType const& value);
                
                /*!
                 * Adds the given choices rewards to this choice.
                 */
                void addRewards(std::vector<ValueType>&& values);
                
                /*!
                 * Retrieves the rewards for this choice.
                 */
                std::vector<ValueType> const& getRewards() const;
                
                /*!
                 * Sets the given values as the rewards of this choice.
                 */
                void setRewards(std::vector<ValueType>&& rewards);
                
                /*!
                 * Resizes the rewards to the desired size and filles newly created values with the provided value.
                 */
                void resizeRewards(std::size_t numberOfRewards, ValueType const& fillValue);
                
                /*!
                 * Retrieves the number of rewards of this choice.
                 */
                std::size_t getNumberOfRewards() const;
                
                /*!
                 * Compresses the underlying distribution.
                 */
                void compress();
                
            private:
                Distribution<IndexType, ValueType>& getMutableDistribution();

                /// The distribution of this choice.
                Distribution<IndexType, ValueType> distribution;
                
                /// The reward values associated with this choice.
                std::vector<ValueType> rewards;
                
                /// A flag storing whether this choice is Markovian.
                bool markovian;
            };
            
        }
    }
}

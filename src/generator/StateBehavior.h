#ifndef STORM_GENERATOR_PRISM_STATEBEHAVIOR_H_
#define STORM_GENERATOR_PRISM_STATEBEHAVIOR_H_

#include <cstdint>

#include "src/generator/Choice.h"

namespace storm {
    namespace generator {
        
        template<typename ValueType, typename StateType = uint32_t>
        class StateBehavior {
        public:
            /*!
             * Creates an empty behavior, i.e. the state was not yet expanded.
             */
            StateBehavior();
            
            /*!
             * Adds the given choice to the behavior of the state.
             */
            void addChoice(Choice<ValueType, StateType>&& choice);
            
            /*!
             * Adds the given state reward to the behavior of the state.
             */
            void addStateReward(ValueType const& stateReward);
            
            /*!
             * Sets whether the state was expanded.
             */
            void setExpanded(bool newValue = true);
            
            /*!
             * Retrieves whether the state was expanded.
             */
            bool wasExpanded() const;
            
            /*!
             * Retrieves whether the behavior is empty in the sense that there are no available choices.
             */
            bool empty() const;
            
            /*!
             * Retrieves an iterator to the choices available in the behavior.
             */
            typename std::vector<Choice<ValueType, StateType>>::const_iterator begin() const;

            /*!
             * Retrieves an iterator past the choices available in the behavior.
             */
            typename std::vector<Choice<ValueType, StateType>>::const_iterator end() const;

            /*!
             * Retrieves the list of state rewards under selected reward models.
             */
            std::vector<ValueType> const& getStateRewards() const;
            
            /*!
             * Retrieves the number of choices in the behavior.
             */
            std::size_t getNumberOfChoices() const;
            
        private:
            // The choices available in the state.
            std::vector<Choice<ValueType, StateType>> choices;
            
            // The state rewards (under the different, selected reward models) of the state.
            std::vector<ValueType> stateRewards;
            
            // A flag indicating whether the state was actually expanded.
            bool expanded;
        };
        
    }
}

#endif /* STORM_GENERATOR_PRISM_STATEBEHAVIOR_H_ */

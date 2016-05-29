#pragma once

#include "src/storage/jani/Variable.h"

namespace storm {
    namespace jani {
        
        class BooleanVariable : public Variable {
        public:
            /*!
             * Creates a boolean variable.
             */
            BooleanVariable(std::string const& name, storm::expressions::Variable const& variable, storm::expressions::Expression const& initialValue = storm::expressions::Expression());
        };
        
    }
}
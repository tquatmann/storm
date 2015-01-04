#include "src/storage/expressions/VariableExpression.h"
#include "src/utility/macros.h"
#include "src/exceptions/InvalidTypeException.h"

namespace storm {
    namespace expressions {
        VariableExpression::VariableExpression(Variable const& variable) : BaseExpression(variable.getManager(), variable.getType()), variable(variable) {
            // Intentionally left empty.
        }
        
        std::string const& VariableExpression::getVariableName() const {
            return variable.getName();
        }
        
        Variable const& VariableExpression::getVariable() const {
            return variable;
        }
        
        bool VariableExpression::evaluateAsBool(Valuation const* valuation) const {
            STORM_LOG_ASSERT(valuation != nullptr, "Evaluating expressions with unknowns without valuation.");
            STORM_LOG_THROW(this->hasBooleanType(), storm::exceptions::InvalidTypeException, "Cannot evaluate expression as boolean: return type is not a boolean.");
            
            return valuation->getBooleanValue(this->getVariable());
        }

        int_fast64_t VariableExpression::evaluateAsInt(Valuation const* valuation) const {
            STORM_LOG_ASSERT(valuation != nullptr, "Evaluating expressions with unknowns without valuation.");
            STORM_LOG_THROW(this->hasIntegralType(), storm::exceptions::InvalidTypeException, "Cannot evaluate expression as integer: return type is not an integer.");
            
            return valuation->getIntegerValue(this->getVariable());
        }
        
        double VariableExpression::evaluateAsDouble(Valuation const* valuation) const {
            STORM_LOG_ASSERT(valuation != nullptr, "Evaluating expressions with unknowns without valuation.");
            STORM_LOG_THROW(this->hasNumericalType(), storm::exceptions::InvalidTypeException, "Cannot evaluate expression as double: return type is not a double.");
            
            if (this->getType().isIntegralType()) {
                return static_cast<double>(valuation->getIntegerValue(this->getVariable()));
            } else {
                return valuation->getRationalValue(this->getVariable());
            }
        }
        
        std::string const& VariableExpression::getIdentifier() const {
            return this->getVariableName();
        }
        
        bool VariableExpression::containsVariables() const {
            return true;
        }

		bool VariableExpression::isVariable() const {
			return true;
		}
        
        std::set<std::string> VariableExpression::getVariables() const {
            return {this->getVariableName()};
        }

        std::shared_ptr<BaseExpression const> VariableExpression::simplify() const {
            return this->shared_from_this();
        }
        
        boost::any VariableExpression::accept(ExpressionVisitor& visitor) const {
            return visitor.visit(*this);
        }
        
        void VariableExpression::printToStream(std::ostream& stream) const {
            stream << this->getVariableName();
        }
    }
}
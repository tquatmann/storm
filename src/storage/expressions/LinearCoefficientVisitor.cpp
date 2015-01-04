#include "src/storage/expressions/LinearCoefficientVisitor.h"

#include "src/storage/expressions/Expressions.h"
#include "src/utility/macros.h"
#include "src/exceptions/InvalidArgumentException.h"

namespace storm {
    namespace expressions {
        LinearCoefficientVisitor::VariableCoefficients::VariableCoefficients(double constantPart) : variableToCoefficientMapping(), constantPart(constantPart) {
            // Intentionally left empty.
        }
        
        LinearCoefficientVisitor::VariableCoefficients& LinearCoefficientVisitor::VariableCoefficients::operator+=(VariableCoefficients&& other) {
            for (auto const& otherVariableCoefficientPair : other.variableToCoefficientMapping) {
                this->variableToCoefficientMapping[otherVariableCoefficientPair.first] += otherVariableCoefficientPair.second;
            }
            constantPart += other.constantPart;
            return *this;
        }
        
        LinearCoefficientVisitor::VariableCoefficients& LinearCoefficientVisitor::VariableCoefficients::operator-=(VariableCoefficients&& other) {
            for (auto const& otherVariableCoefficientPair : other.variableToCoefficientMapping) {
                this->variableToCoefficientMapping[otherVariableCoefficientPair.first] -= otherVariableCoefficientPair.second;
            }
            constantPart -= other.constantPart;
            return *this;
        }
        
        LinearCoefficientVisitor::VariableCoefficients& LinearCoefficientVisitor::VariableCoefficients::operator*=(VariableCoefficients&& other) {
            STORM_LOG_THROW(variableToCoefficientMapping.size() == 0 || other.variableToCoefficientMapping.size() == 0, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
            if (other.variableToCoefficientMapping.size() > 0) {
                variableToCoefficientMapping = std::move(other.variableToCoefficientMapping);
                std::swap(constantPart, other.constantPart);
            }
            for (auto const& otherVariableCoefficientPair : other.variableToCoefficientMapping) {
                this->variableToCoefficientMapping[otherVariableCoefficientPair.first] *= other.constantPart;
            }
            constantPart *= other.constantPart;
            return *this;
        }

        LinearCoefficientVisitor::VariableCoefficients& LinearCoefficientVisitor::VariableCoefficients::operator/=(VariableCoefficients&& other) {
            STORM_LOG_THROW(other.variableToCoefficientMapping.size() == 0, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
            for (auto const& otherVariableCoefficientPair : other.variableToCoefficientMapping) {
                this->variableToCoefficientMapping[otherVariableCoefficientPair.first] /= other.constantPart;
            }
            constantPart /= other.constantPart;
            return *this;
        }
        
        void LinearCoefficientVisitor::VariableCoefficients::negate() {
            for (auto& variableCoefficientPair : variableToCoefficientMapping) {
                variableCoefficientPair.second = -variableCoefficientPair.second;
            }
            constantPart = -constantPart;
        }

        void LinearCoefficientVisitor::VariableCoefficients::setCoefficient(storm::expressions::Variable const& variable, double coefficient) {
            variableToCoefficientMapping[variable] = coefficient;
        }
        
        double LinearCoefficientVisitor::VariableCoefficients::getCoefficient(storm::expressions::Variable const& variable) {
            return variableToCoefficientMapping[variable];
        }
        
        LinearCoefficientVisitor::VariableCoefficients LinearCoefficientVisitor::getLinearCoefficients(Expression const& expression) {
            return boost::any_cast<VariableCoefficients>(expression.getBaseExpression().accept(*this));
        }
        
        boost::any LinearCoefficientVisitor::visit(IfThenElseExpression const& expression) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
        }
        
        boost::any LinearCoefficientVisitor::visit(BinaryBooleanFunctionExpression const& expression) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
        }
        
        boost::any LinearCoefficientVisitor::visit(BinaryNumericalFunctionExpression const& expression) {
            VariableCoefficients leftResult = boost::any_cast<VariableCoefficients>(expression.getFirstOperand()->accept(*this));
            VariableCoefficients rightResult = boost::any_cast<VariableCoefficients>(expression.getSecondOperand()->accept(*this));

            if (expression.getOperatorType() == BinaryNumericalFunctionExpression::OperatorType::Plus) {
                leftResult += std::move(rightResult);
            } else if (expression.getOperatorType() == BinaryNumericalFunctionExpression::OperatorType::Minus) {
                leftResult -= std::move(rightResult);
            } else if (expression.getOperatorType() == BinaryNumericalFunctionExpression::OperatorType::Times) {
                leftResult *= std::move(rightResult);
            } else if (expression.getOperatorType() == BinaryNumericalFunctionExpression::OperatorType::Divide) {
                leftResult /= std::move(rightResult);
            } else {
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
            }
            return rightResult;
        }
        
        boost::any LinearCoefficientVisitor::visit(BinaryRelationExpression const& expression) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
        }
        
        boost::any LinearCoefficientVisitor::visit(VariableExpression const& expression) {
            VariableCoefficients coefficients;
            if (expression.getType().isNumericalType()) {
                coefficients.setCoefficient(expression.getVariable(), 1);
            } else {
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
            }
            return coefficients;
        }
        
        boost::any LinearCoefficientVisitor::visit(UnaryBooleanFunctionExpression const& expression) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
        }
        
        boost::any LinearCoefficientVisitor::visit(UnaryNumericalFunctionExpression const& expression) {
            VariableCoefficients childResult = boost::any_cast<VariableCoefficients>(expression.getOperand()->accept(*this));
            
            if (expression.getOperatorType() == UnaryNumericalFunctionExpression::OperatorType::Minus) {
                childResult.negate();
                return childResult;
            } else {
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
            }
        }
        
        boost::any LinearCoefficientVisitor::visit(BooleanLiteralExpression const& expression) {
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Expression is non-linear.");
        }
        
        boost::any LinearCoefficientVisitor::visit(IntegerLiteralExpression const& expression) {
            return VariableCoefficients(static_cast<double>(expression.getValue()));
        }
        
        boost::any LinearCoefficientVisitor::visit(DoubleLiteralExpression const& expression) {
            return VariableCoefficients(expression.getValue());
        }
    }
}
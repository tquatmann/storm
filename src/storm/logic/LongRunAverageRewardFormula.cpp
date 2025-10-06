#include "storm/logic/LongRunAverageRewardFormula.h"
#include <boost/any.hpp>
#include <ostream>

#include "storm/logic/FormulaVisitor.h"

namespace storm {
namespace logic {
LongRunAverageRewardFormula::LongRunAverageRewardFormula(boost::optional<RewardAccumulation> rewardAccumulation) : rewardAccumulation(rewardAccumulation) {
    // Intentionally left empty.
}

LongRunAverageRewardFormula::LongRunAverageRewardFormula(std::string const& boundRewardModelName, storm::logic::Bound const& bound,
                                                         boost::optional<RewardAccumulation> rewardAccumulation)
    : PathFormula(), rewardAccumulation(rewardAccumulation), boundRewardModelName(boundRewardModelName), bound(bound) {
    // Intentionally left empty.
}

bool LongRunAverageRewardFormula::isLongRunAverageRewardFormula() const {
    return true;
}

bool LongRunAverageRewardFormula::isRewardPathFormula() const {
    return !hasBound();
}

bool LongRunAverageRewardFormula::isProbabilityPathFormula() const {
    return hasBound();
}

bool LongRunAverageRewardFormula::hasRewardAccumulation() const {
    return rewardAccumulation.is_initialized();
}

RewardAccumulation const& LongRunAverageRewardFormula::getRewardAccumulation() const {
    return rewardAccumulation.get();
}

std::shared_ptr<LongRunAverageRewardFormula const> LongRunAverageRewardFormula::stripRewardAccumulation() const {
    if (hasBound()) {
        return std::make_shared<LongRunAverageRewardFormula>(getBoundRewardModelName(), getBound());

    } else {
        return std::make_shared<LongRunAverageRewardFormula>();
    }
}

bool LongRunAverageRewardFormula::hasBound() const {
    STORM_LOG_ASSERT(bound.has_value() == boundRewardModelName.has_value(), "inconsistent state");
    return bound.has_value();
}

storm::logic::Bound const& LongRunAverageRewardFormula::getBound() const {
    STORM_LOG_ASSERT(hasBound(), "Tried to retrieve a bound but this long-run average reward formula does not have a bound.");
    return bound.value();
}

void LongRunAverageRewardFormula::gatherReferencedRewardModels(std::set<std::string>& referencedRewardModels) const {
    if (this->hasBound()) {
        referencedRewardModels.insert(this->getBoundRewardModelName());
    }
}

void LongRunAverageRewardFormula::gatherUsedVariables(std::set<storm::expressions::Variable>& usedVariables) const {
    PathFormula::gatherUsedVariables(usedVariables);
    if (this->hasBound()) {
        getBound().threshold.gatherVariables(usedVariables);
    }
}

std::string const& LongRunAverageRewardFormula::getBoundRewardModelName() const {
    return boundRewardModelName.value();
}

boost::any LongRunAverageRewardFormula::accept(FormulaVisitor const& visitor, boost::any const& data) const {
    return visitor.visit(*this, data);
}

std::ostream& LongRunAverageRewardFormula::writeToStream(std::ostream& out, bool /* allowParentheses */) const {
    // No parentheses necessary
    out << (hasBound() ? "LRASAT" : "LRA");
    if (hasRewardAccumulation()) {
        out << "[" << getRewardAccumulation() << "]";
    }
    if (hasBound()) {
        out << "{\"" << getBoundRewardModelName() << "\"}" << getBound().comparisonType << getBound().threshold;
    }
    return out;
}

}  // namespace logic
}  // namespace storm

#ifndef STORM_LOGIC_LONGRUNAVERAGEREWARDFORMULA_H_
#define STORM_LOGIC_LONGRUNAVERAGEREWARDFORMULA_H_

#include <optional>
#include <string>
#include "storm/logic/Bound.h"
#include "storm/logic/PathFormula.h"
#include "storm/logic/RewardAccumulation.h"

namespace storm {
namespace logic {
class LongRunAverageRewardFormula : public PathFormula {
   public:
    LongRunAverageRewardFormula(boost::optional<RewardAccumulation> rewardAccumulation = boost::none);
    LongRunAverageRewardFormula(std::string const& boundRewardModelName, storm::logic::Bound const& bound,
                                boost::optional<RewardAccumulation> rewardAccumulation = boost::none);

    virtual ~LongRunAverageRewardFormula() {
        // Intentionally left empty.
    }

    virtual bool isLongRunAverageRewardFormula() const override;
    virtual bool isRewardPathFormula() const override;
    virtual bool isProbabilityPathFormula() const override;
    bool hasRewardAccumulation() const;
    RewardAccumulation const& getRewardAccumulation() const;
    std::shared_ptr<LongRunAverageRewardFormula const> stripRewardAccumulation() const;
    bool hasBound() const;
    storm::logic::Bound const& getBound() const;
    std::string const& getBoundRewardModelName() const;

    virtual void gatherReferencedRewardModels(std::set<std::string>& referencedRewardModels) const override;
    virtual void gatherUsedVariables(std::set<storm::expressions::Variable>& usedVariables) const override;

    virtual boost::any accept(FormulaVisitor const& visitor, boost::any const& data) const override;

    virtual std::ostream& writeToStream(std::ostream& out, bool allowParentheses = false) const override;

   private:
    boost::optional<RewardAccumulation> rewardAccumulation;
    std::optional<std::string> boundRewardModelName;
    std::optional<storm::logic::Bound> bound;
};
}  // namespace logic
}  // namespace storm

#endif /* STORM_LOGIC_LONGRUNAVERAGEREWARDFORMULA_H_ */

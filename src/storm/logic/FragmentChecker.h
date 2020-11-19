#ifndef STORM_LOGIC_FRAGMENTCHECKER_H_
#define STORM_LOGIC_FRAGMENTCHECKER_H_

#include "storm/logic/FormulaVisitor.h"

#include "storm/logic/FragmentSpecification.h"

namespace storm {
    namespace logic {
        
        class FragmentChecker : public FormulaVisitor {
        public:
            bool conformsToSpecification(Formula const& f, FragmentSpecification const& specification) const;
            
            virtual boost::any visit(AtomicExpressionFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(AtomicLabelFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(BinaryBooleanStateFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(BooleanLiteralFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(BoundedUntilFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(ConditionalFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(CumulativeRewardFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(EventuallyFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(TimeOperatorFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(GloballyFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(GameFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(InstantaneousRewardFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(LongRunAverageOperatorFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(LongRunAverageRewardFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(MultiObjectiveFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(QuantileFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(NextFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(ProbabilityOperatorFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(RewardOperatorFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(TotalRewardFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(UnaryBooleanStateFormula const& f, boost::any const& data) const override;
            virtual boost::any visit(UntilFormula const& f, boost::any const& data) const override;
        };
        
    }
}


#endif /* STORM_LOGIC_FRAGMENTCHECKER_H_ */

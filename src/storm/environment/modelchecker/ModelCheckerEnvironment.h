#pragma once

#include <boost/optional.hpp>
#include <memory>
#include <string>

#include "storm/automata/AutomatonType.h"
#include "storm/environment/Environment.h"
#include "storm/environment/SubEnvironment.h"
#include "storm/modelchecker/helper/conditional/ConditionalAlgorithmSetting.h"
#include "storm/modelchecker/helper/infinitehorizon/SteadyStateDistributionAlgorithm.h"

namespace storm {

// Forward declare subenvironments
class MultiObjectiveModelCheckerEnvironment;

class ModelCheckerEnvironment {
   public:
    ModelCheckerEnvironment();
    ~ModelCheckerEnvironment();

    MultiObjectiveModelCheckerEnvironment& multi();
    MultiObjectiveModelCheckerEnvironment const& multi() const;

    SteadyStateDistributionAlgorithm getSteadyStateDistributionAlgorithm() const;
    void setSteadyStateDistributionAlgorithm(SteadyStateDistributionAlgorithm value);

    ConditionalAlgorithmSetting getConditionalAlgorithmSetting() const;
    void setConditionalAlgorithmSetting(ConditionalAlgorithmSetting value);

    bool isLtl2AutToolSet() const;
    std::string const& getLtl2AutTool() const;
    void setLtl2AutTool(std::string const& value);
    void unsetLtl2AutTool();
    storm::automata::AutomatonType getLtlAutomatonType() const;
    void setLtlAutomatonType(storm::automata::AutomatonType const& type);

   private:
    SubEnvironment<MultiObjectiveModelCheckerEnvironment> multiObjectiveModelCheckerEnvironment;
    std::optional<std::string> ltl2AutTool;
    storm::automata::AutomatonType ltlAutomatonType;
    SteadyStateDistributionAlgorithm steadyStateDistributionAlgorithm;
    ConditionalAlgorithmSetting conditionalAlgorithmSetting;
};
}  // namespace storm

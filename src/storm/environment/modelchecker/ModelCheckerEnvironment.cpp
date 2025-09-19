#include "storm/environment/modelchecker/ModelCheckerEnvironment.h"

#include "storm/environment/modelchecker/MultiObjectiveModelCheckerEnvironment.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/IOSettings.h"
#include "storm/settings/modules/ModelCheckerSettings.h"
#include "storm/utility/macros.h"

#include "storm/exceptions/InvalidEnvironmentException.h"
#include "storm/exceptions/UnexpectedException.h"

namespace storm {

ModelCheckerEnvironment::ModelCheckerEnvironment() {
    auto const& mcSettings = storm::settings::getModule<storm::settings::modules::ModelCheckerSettings>();
    if (mcSettings.isLtl2AutToolSet()) {
        ltl2AutTool = mcSettings.getLtl2AutTool();
        ltlAutomatonType = mcSettings.getLtlAutomatonType();
    }
    auto const& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    steadyStateDistributionAlgorithm = ioSettings.getSteadyStateDistributionAlgorithm();

    conditionalAlgorithmSetting = mcSettings.getConditionalAlgorithmSetting();
}

ModelCheckerEnvironment::~ModelCheckerEnvironment() {
    // Intentionally left empty
}

SteadyStateDistributionAlgorithm ModelCheckerEnvironment::getSteadyStateDistributionAlgorithm() const {
    return steadyStateDistributionAlgorithm;
}

void ModelCheckerEnvironment::setSteadyStateDistributionAlgorithm(SteadyStateDistributionAlgorithm value) {
    steadyStateDistributionAlgorithm = value;
}

ConditionalAlgorithmSetting ModelCheckerEnvironment::getConditionalAlgorithmSetting() const {
    return conditionalAlgorithmSetting;
}

void ModelCheckerEnvironment::setConditionalAlgorithmSetting(ConditionalAlgorithmSetting value) {
    conditionalAlgorithmSetting = value;
}

MultiObjectiveModelCheckerEnvironment& ModelCheckerEnvironment::multi() {
    return multiObjectiveModelCheckerEnvironment.get();
}

MultiObjectiveModelCheckerEnvironment const& ModelCheckerEnvironment::multi() const {
    return multiObjectiveModelCheckerEnvironment.get();
}

bool ModelCheckerEnvironment::isLtl2AutToolSet() const {
    return ltl2AutTool.has_value();
}

std::string const& ModelCheckerEnvironment::getLtl2AutTool() const {
    return ltl2AutTool.value();
}

void ModelCheckerEnvironment::setLtl2AutTool(std::string const& value) {
    ltl2AutTool = value;
}

void ModelCheckerEnvironment::unsetLtl2AutTool() {
    ltl2AutTool = std::nullopt;
}

storm::automata::AutomatonType ModelCheckerEnvironment::getLtlAutomatonType() const {
    return ltlAutomatonType;
}

void ModelCheckerEnvironment::setLtlAutomatonType(storm::automata::AutomatonType const& type) {
    ltlAutomatonType = type;
}
}  // namespace storm

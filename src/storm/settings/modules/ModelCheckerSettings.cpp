#include "storm/settings/modules/ModelCheckerSettings.h"

#include "storm/settings/Argument.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingsManager.h"

namespace storm::settings::modules {

const std::string ModelCheckerSettings::moduleName = "modelchecker";
const std::string ModelCheckerSettings::filterRewZeroOptionName = "filterrewzero";
const std::string ModelCheckerSettings::ltl2AutToolOptionName = "ltl2auttool";
const std::string ModelCheckerSettings::ltlAutomatonType = "ltlauttype";
const std::string ModelCheckerSettings::conditionalAlgorithmOptionName = "conditional";

ModelCheckerSettings::ModelCheckerSettings() : ModuleSettings(moduleName) {
    this->addOption(storm::settings::OptionBuilder(moduleName, filterRewZeroOptionName, false,
                                                   "If set, states with reward zero are filtered out, potentially reducing the size of the equation system")
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, ltl2AutToolOptionName, false,
                                                   "If set, use an external tool to convert LTL formulas to state-based automata in HOA format")
                        .setIsAdvanced()
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument(
                                         "filename", "A script that can be called with a formula and a name for the output automaton.")
                                         .build())
                        .build());

    std::vector<std::string> const automatontypes = {"da", "ldba"};
    this->addOption(storm::settings::OptionBuilder(moduleName, ltlAutomatonType, false, "The type of automaton used for ltl model checking.")
                        .setIsAdvanced()
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("type", "The desired type of automaton.")
                                         .addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(automatontypes))
                                         .setDefaultValueString("da")
                                         .build())
                        .build());

    std::vector<std::string> const conditionalAlgs = {"default", "restart", "bisection", "bisection-advanced", "pi"};
    this->addOption(storm::settings::OptionBuilder(moduleName, conditionalAlgorithmOptionName, false, "The used algorithm for conditional probabilities.")
                        .setIsAdvanced()
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the method to use.")
                                         .addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(conditionalAlgs))
                                         .setDefaultValueString("default")
                                         .build())
                        .build());
}

bool ModelCheckerSettings::isFilterRewZeroSet() const {
    return this->getOption(filterRewZeroOptionName).getHasOptionBeenSet();
}

bool ModelCheckerSettings::isLtl2AutToolSet() const {
    return this->getOption(ltl2AutToolOptionName).getHasOptionBeenSet();
}

std::string ModelCheckerSettings::getLtl2AutTool() const {
    return this->getOption(ltl2AutToolOptionName).getArgumentByName("filename").getValueAsString();
}

storm::automata::AutomatonType ModelCheckerSettings::getLtlAutomatonType() const {
    return storm::automata::automatonTypeFromString(this->getOption(ltlAutomatonType).getArgumentByName("type").getValueAsString());
}

bool ModelCheckerSettings::isConditionalAlgorithmSet() const {
    return this->getOption(conditionalAlgorithmOptionName).getHasOptionBeenSet();
}

ConditionalAlgorithmSetting ModelCheckerSettings::getConditionalAlgorithmSetting() const {
    return conditionalAlgorithmSettingFromString(this->getOption(conditionalAlgorithmOptionName).getArgumentByName("name").getValueAsString());
}

}  // namespace storm::settings::modules

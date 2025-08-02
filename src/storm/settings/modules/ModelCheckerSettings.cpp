#include "storm/settings/modules/ModelCheckerSettings.h"

#include "storm/settings/Argument.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/SettingsManager.h"

namespace storm {
namespace settings {
namespace modules {

const std::string ModelCheckerSettings::moduleName = "modelchecker";
const std::string ModelCheckerSettings::filterRewZeroOptionName = "filterrewzero";
const std::string ModelCheckerSettings::ltl2daToolOptionName = "ltl2datool";
const std::string ModelCheckerSettings::ltl2daSyntax = "ltl2dasyntax";

ModelCheckerSettings::ModelCheckerSettings() : ModuleSettings(moduleName) {
    this->addOption(storm::settings::OptionBuilder(moduleName, filterRewZeroOptionName, false,
                                                   "If set, states with reward zero are filtered out, potentially reducing the size of the equation system")
                        .setIsAdvanced()
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, ltl2daToolOptionName, false,
                                                   "If set, use an external tool to convert LTL formulas to state-based deterministic automata in HOA format")
                        .setIsAdvanced()
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument(
                                         "filename", "A script that can be called with a prefix formula and a name for the output automaton.")
                                         .build())
                        .build());
    this->addOption(storm::settings::OptionBuilder(moduleName, ltl2daSyntax, false,
                                                   "Provides information about the type of automaton used, either deterministic or limit-deterministic.")
                        .setIsAdvanced()
                        .addArgument(storm::settings::ArgumentBuilder::createStringArgument(
                                        "syntax", "Either 'deterministic' or 'limit-deterministic', depending on the desired type of automaton.")
                                        .build())
                        .build());
}

bool ModelCheckerSettings::isFilterRewZeroSet() const {
    return this->getOption(filterRewZeroOptionName).getHasOptionBeenSet();
}

bool ModelCheckerSettings::isLtl2daToolSet() const {
    return this->getOption(ltl2daToolOptionName).getHasOptionBeenSet();
}

std::string ModelCheckerSettings::getLtl2daTool() const {
    return this->getOption(ltl2daToolOptionName).getArgumentByName("filename").getValueAsString();
}

bool ModelCheckerSettings::isLtl2DeterministicAutomatonSet() const {
    return this->getOption(ltl2daSyntax).getArgumentByName("syntax").getValueAsString() == "deterministic";
}


}  // namespace modules
}  // namespace settings
}  // namespace storm

#include "storm/settings/modules/LearningSettings.h"

#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/SettingsManager.h"

namespace storm {
namespace settings {
namespace modules {

const std::string LearningSettings::moduleName = "learning";
const std::string LearningSettings::learnIMDPFromMDPOptionName = "learn-imdp-from-mdp";

LearningSettings::LearningSettings() : ModuleSettings(moduleName) {
    this->addOption(
        storm::settings::OptionBuilder(moduleName, learnIMDPFromMDPOptionName, false,
                                       "Learns an IMDP by computing Clopper-Pearson intervals based on the input MDP.")
            .setIsAdvanced()
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("delta", "Specifies the confidence interval by 1 - delta.")
                             .setDefaultValueDouble(0.01)
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createIntegerArgument("number-of-samples", "Specifies the number of samples per state-action pair.")
                             .setDefaultValueInteger(1000000)
                             .addValidatorInteger(ArgumentValidatorFactory::createIntegerGreaterValidator(0))
                             .build())
            .build());
}

bool LearningSettings::isLearnIMDPFromMDPSet() {
    return this->getOption(learnIMDPFromMDPOptionName).getHasOptionBeenSet();
}

double LearningSettings::getDeltaValue() {
    return this->getOption(learnIMDPFromMDPOptionName).getArgumentByName("delta").getValueAsDouble();
}

uint_fast64_t LearningSettings::getNumberOfSamples() {
    return this->getOption(learnIMDPFromMDPOptionName).getArgumentByName("number-of-samples").getValueAsInteger();
}

void LearningSettings::finalize() {
    // Intentionally left empty.
}

bool LearningSettings::check() const {
    return true;
}

}  // namespace modules
}  // namespace settings
}  // namespace storm
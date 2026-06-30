#include "storm/settings/modules/LearningSettings.h"

#include "storm/settings/ArgumentBuilder.h"
#include "storm/exceptions/InvalidSettingsException.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/SettingsManager.h"

namespace storm {
namespace settings {
namespace modules {

const std::string LearningSettings::moduleName = "learning";
const std::string LearningSettings::learnIMDPOptionName = "learn-imdp";
const std::string LearningSettings::learnIMDPUntilWidthOptionName = "learn-imdp-until-width";
const std::string LearningSettings::learnIMDPUntilL1WidthOptionName = "learn-imdp-until-l1-width";

LearningSettings::LearningSettings() : ModuleSettings(moduleName) {
    this->addOption(
        storm::settings::OptionBuilder(moduleName, learnIMDPOptionName, false,
                                       "Learns an IMDP by computing Clopper-Pearson intervals based on the input model.")
            .setIsAdvanced()
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("delta", "Specifies the confidence interval by 1 - delta.")
                             .setDefaultValueDouble(0.01)
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createIntegerArgument("number-of-samples", "Specifies the number of samples per state-action pair.")
                             .setDefaultValueInteger(1000000)
                             .addValidatorInteger(ArgumentValidatorFactory::createIntegerGreaterValidator(0))
                             .build())
            .build());

    this->addOption(
        storm::settings::OptionBuilder(moduleName, learnIMDPUntilWidthOptionName, false,
                                       "Learns an IMDP by sampling until all Clopper-Pearson intervals are below a given width.")
            .setIsAdvanced()
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("delta", "Specifies the confidence interval by 1 - delta.")
                             .setDefaultValueDouble(0.01)
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("max-width", "Specifies the maximum allowed interval width.")
                             .addValidatorDouble(ArgumentValidatorFactory::createDoubleGreaterValidator(0.0))
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createIntegerArgument("max-number-of-samples", "Specifies the maximum number of samples per state-action pair.")
                             .setDefaultValueInteger(1000000)
                             .addValidatorInteger(ArgumentValidatorFactory::createIntegerGreaterValidator(0))
                             .build())
            .build());

    this->addOption(
        storm::settings::OptionBuilder(moduleName, learnIMDPUntilL1WidthOptionName, false,
                                       "Learns an IMDP by sampling until the feasible L1 diameter of each state-action successor confidence polytope is below a given threshold.")
            .setIsAdvanced()
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("delta", "Specifies the confidence interval by 1 - delta.")
                             .setDefaultValueDouble(0.01)
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("max-l1-width", "Specifies the maximum allowed feasible L1 diameter per state-action successor distribution.")
                             .addValidatorDouble(ArgumentValidatorFactory::createDoubleGreaterValidator(0.0))
                             .build())
            .addArgument(storm::settings::ArgumentBuilder::createIntegerArgument("max-number-of-samples", "Specifies the maximum number of samples per state-action pair.")
                             .setDefaultValueInteger(1000000)
                             .addValidatorInteger(ArgumentValidatorFactory::createIntegerGreaterValidator(0))
                             .build())
            .build());
}

bool LearningSettings::isLearnIMDPFromMDPSet() {
    return this->getOption(learnIMDPOptionName).getHasOptionBeenSet();
}

bool LearningSettings::isLearnIMDPFromMDPWithMaxWidthSet() {
    return this->getOption(learnIMDPUntilWidthOptionName).getHasOptionBeenSet();
}

bool LearningSettings::isLearnIMDPFromMDPWithMaxL1WidthSet() {
    return this->getOption(learnIMDPUntilL1WidthOptionName).getHasOptionBeenSet();
}

double LearningSettings::getDeltaValue() {
    return this->getOption(learnIMDPOptionName).getArgumentByName("delta").getValueAsDouble();
}

uint_fast64_t LearningSettings::getNumberOfSamples() {
    return this->getOption(learnIMDPOptionName).getArgumentByName("number-of-samples").getValueAsInteger();
}

double LearningSettings::getWidthDeltaValue() {
    return this->getOption(learnIMDPUntilWidthOptionName).getArgumentByName("delta").getValueAsDouble();
}

double LearningSettings::getMaxWidth() {
    return this->getOption(learnIMDPUntilWidthOptionName).getArgumentByName("max-width").getValueAsDouble();
}

uint_fast64_t LearningSettings::getMaxNumberOfSamples() {
    return this->getOption(learnIMDPUntilWidthOptionName).getArgumentByName("max-number-of-samples").getValueAsInteger();
}

double LearningSettings::getL1WidthDeltaValue() {
    return this->getOption(learnIMDPUntilL1WidthOptionName).getArgumentByName("delta").getValueAsDouble();
}

double LearningSettings::getMaxL1Width() {
    return this->getOption(learnIMDPUntilL1WidthOptionName).getArgumentByName("max-l1-width").getValueAsDouble();
}

uint_fast64_t LearningSettings::getL1WidthMaxNumberOfSamples() {
    return this->getOption(learnIMDPUntilL1WidthOptionName).getArgumentByName("max-number-of-samples").getValueAsInteger();
}

void LearningSettings::finalize() {
    // Intentionally left empty.
}

bool LearningSettings::check() const {
    auto numberOfLearningOptionsSet = static_cast<uint_fast64_t>(this->getOption(learnIMDPOptionName).getHasOptionBeenSet()) +
                                      static_cast<uint_fast64_t>(this->getOption(learnIMDPUntilWidthOptionName).getHasOptionBeenSet()) +
                                      static_cast<uint_fast64_t>(this->getOption(learnIMDPUntilL1WidthOptionName).getHasOptionBeenSet());
    STORM_LOG_THROW(numberOfLearningOptionsSet <= 1, storm::exceptions::InvalidSettingsException,
                    "Cannot set more than one learning option.");
    return true;
}

}  // namespace modules
}  // namespace settings
}  // namespace storm

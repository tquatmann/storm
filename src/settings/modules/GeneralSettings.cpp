#include "src/settings/modules/GeneralSettings.h"

#include "src/settings/SettingsManager.h"
#include "src/settings/SettingMemento.h"
#include "src/settings/Option.h"
#include "src/settings/OptionBuilder.h"
#include "src/settings/ArgumentBuilder.h"
#include "src/settings/Argument.h"
#include "src/solver/SolverSelectionOptions.h"

#include "src/storage/dd/DdType.h"

#include "src/exceptions/InvalidSettingsException.h"

namespace storm {
    namespace settings {
        namespace modules {
            
            const std::string GeneralSettings::moduleName = "general";
            const std::string GeneralSettings::helpOptionName = "help";
            const std::string GeneralSettings::helpOptionShortName = "h";
            const std::string GeneralSettings::printTimingsOptionName = "printTime";
            const std::string GeneralSettings::printTimingsOptionShortName = "pt";
            const std::string GeneralSettings::versionOptionName = "version";
            const std::string GeneralSettings::verboseOptionName = "verbose";
            const std::string GeneralSettings::verboseOptionShortName = "v";
            const std::string GeneralSettings::precisionOptionName = "precision";
            const std::string GeneralSettings::precisionOptionShortName = "eps";
            const std::string GeneralSettings::configOptionName = "config";
            const std::string GeneralSettings::configOptionShortName = "c";
            const std::string GeneralSettings::propertyOptionName = "prop";
            const std::string GeneralSettings::propertyOptionShortName = "prop";
            const std::string GeneralSettings::bisimulationOptionName = "bisimulation";
            const std::string GeneralSettings::bisimulationOptionShortName = "bisim";
            const std::string GeneralSettings::parametricOptionName = "parametric";
            const std::string GeneralSettings::parametricRegionOptionName = "parametricRegion";
            const std::string GeneralSettings::exactOptionName = "exact";


            GeneralSettings::GeneralSettings() : ModuleSettings(moduleName) {
                this->addOption(storm::settings::OptionBuilder(moduleName, helpOptionName, false, "Shows all available options, arguments and descriptions.").setShortName(helpOptionShortName)
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("hint", "A regular expression to show help for all matching entities or 'all' for the complete help.").setDefaultValueString("all").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, versionOptionName, false, "Prints the version information.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, printTimingsOptionName, false, "Prints the timing at the end").setShortName(printTimingsOptionShortName).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, verboseOptionName, false, "Enables more verbose output.").setShortName(verboseOptionShortName).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, precisionOptionName, false, "The internally used precision.").setShortName(precisionOptionShortName)
                                .addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("value", "The precision to use.").setDefaultValueDouble(1e-06).addValidationFunctionDouble(storm::settings::ArgumentValidators::doubleRangeValidatorExcluding(0.0, 1.0)).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, configOptionName, false, "If given, this file will be read and parsed for additional configuration settings.").setShortName(configOptionShortName)
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("filename", "The name of the file from which to read the configuration.").addValidationFunctionString(storm::settings::ArgumentValidators::existingReadableFileValidator()).build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, propertyOptionName, false, "Specifies the formulas to be checked on the model.").setShortName(propertyOptionShortName)
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("formula or filename", "The formula or the file containing the formulas.").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, parametricRegionOptionName, false, "Sets whether to use the parametric Region engine.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, bisimulationOptionName, false, "Sets whether to perform bisimulation minimization.").setShortName(bisimulationOptionShortName).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, parametricOptionName, false, "Sets whether to enable parametric model checking.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exactOptionName, false, "Sets whether to enable exact model checking.").build());
            }
            
            bool GeneralSettings::isHelpSet() const {
                return this->getOption(helpOptionName).getHasOptionBeenSet();
            }
            
            bool GeneralSettings::isVersionSet() const {
                return this->getOption(versionOptionName).getHasOptionBeenSet();
            }
            
            std::string GeneralSettings::getHelpModuleName() const {
                return this->getOption(helpOptionName).getArgumentByName("hint").getValueAsString();
            }
            
            bool GeneralSettings::isVerboseSet() const {
                return this->getOption(verboseOptionName).getHasOptionBeenSet();
            }
            
            double GeneralSettings::getPrecision() const {
                return this->getOption(precisionOptionName).getArgumentByName("value").getValueAsDouble();
            }
            
            bool GeneralSettings::isConfigSet() const {
                return this->getOption(configOptionName).getHasOptionBeenSet();
            }
            
            std::string GeneralSettings::getConfigFilename() const {
                return this->getOption(configOptionName).getArgumentByName("filename").getValueAsString();
            }
            
            bool GeneralSettings::isPropertySet() const {
                return this->getOption(propertyOptionName).getHasOptionBeenSet();
            }
            
            std::string GeneralSettings::getProperty() const {
                return this->getOption(propertyOptionName).getArgumentByName("formula or filename").getValueAsString();
            }
                                    
            bool GeneralSettings::isBisimulationSet() const {
                return this->getOption(bisimulationOptionName).getHasOptionBeenSet();
            }
            
            bool GeneralSettings::isParametricSet() const {
                return this->getOption(parametricOptionName).getHasOptionBeenSet();
            }

            bool GeneralSettings::isParametricRegionSet() const {
                return this->getOption(parametricRegionOptionName).getHasOptionBeenSet();
            }
                                
            bool GeneralSettings::isPrintTimingsSet() const {
                return this->getOption(printTimingsOptionName).getHasOptionBeenSet();
            }

            bool GeneralSettings::isExactSet() const {
                return this->getOption(exactOptionName).getHasOptionBeenSet();
            }
            
            void GeneralSettings::finalize() {
            }

            bool GeneralSettings::check() const {
                return true;
            }
            
        } // namespace modules
    } // namespace settings
} // namespace storm

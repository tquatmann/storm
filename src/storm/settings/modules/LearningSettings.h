#pragma once

#include "storm/settings/modules/ModuleSettings.h"

namespace storm {
namespace settings {
namespace modules {

class LearningSettings : public ModuleSettings {
   public:
    LearningSettings();

    bool isLearnIMDPFromMDPSet();
    bool isLearnIMDPFromMDPWithMaxL1WidthSet();

    double getLambdaValue();

    uint_fast64_t getNumberOfSamples();

    double getL1WidthLambdaValue();
    double getMaxL1Width();
    uint_fast64_t getL1WidthMaxNumberOfSamples();

    bool check() const override;
    void finalize() override;

    // The name of the module.
    static const std::string moduleName;

   private:
    // Define the string names of the options as constants.
    static const std::string learnIMDPOptionName;
    static const std::string learnIMDPUntilL1WidthOptionName;
};

}  // namespace modules
}  // namespace settings
}  // namespace storm

#pragma once

#include "storm/settings/modules/ModuleSettings.h"

namespace storm {
namespace settings {
namespace modules {

class LearningSettings : public ModuleSettings {
   public:
    LearningSettings();

    bool isLearnIMDPFromMDPSet();

    double getDeltaValue();

    uint_fast64_t getNumberOfSamples();

    bool check() const override;
    void finalize() override;

    // The name of the module.
    static const std::string moduleName;

   private:
    // Define the string names of the options as constants.
    static const std::string learnIMDPFromMDPOptionName;
};

}  // namespace modules
}  // namespace settings
}  // namespace storm
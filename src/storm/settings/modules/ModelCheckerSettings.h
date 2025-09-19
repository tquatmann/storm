#pragma once

#include "storm/automata/AutomatonType.h"
#include "storm/modelchecker/helper/conditional/ConditionalAlgorithmSetting.h"
#include "storm/settings/modules/ModuleSettings.h"

namespace storm {
namespace settings {
namespace modules {

/*!
 * This class represents the general settings.
 */
class ModelCheckerSettings : public ModuleSettings {
   public:
    /*!
     * Creates a new set of general settings.
     */
    ModelCheckerSettings();

    bool isFilterRewZeroSet() const;

    /*!
     * Retrieves whether the external ltl2Aut tool has been set.
     *
     * @return True iff the external ltl2Aut tool has been set.
     */
    bool isLtl2AutToolSet() const;

    /*!
     * Retrieves the external ltl2aut tool that is used for converting LTL formulas to automata.
     *
     * @return The executable to use for converting LTL formulas to automata.
     */
    std::string getLtl2AutTool() const;

    /*!
     * Retrieves whether the deterministic or limit-deterministic automata are used for the LTL model checking.
     *
     * @return True if deterministic Rabin automata are to be used
     */
    storm::automata::AutomatonType getLtlAutomatonType() const;

    /*!
     * Retrieves whether an algorithm for conditional properties has been set.
     */
    bool isConditionalAlgorithmSet() const;

    /*!
     * Retrieves the specified algorithm for conditional probabilities.
     */
    ConditionalAlgorithmSetting getConditionalAlgorithmSetting() const;

    // The name of the module.
    static const std::string moduleName;

   private:
    // Define the string names of the options as constants.
    static const std::string filterRewZeroOptionName;
    static const std::string ltl2AutToolOptionName;
    static const std::string ltlAutomatonType;
    static const std::string conditionalAlgorithmOptionName;
};

}  // namespace modules
}  // namespace settings
}  // namespace storm

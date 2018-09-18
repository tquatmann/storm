#pragma once

#include <string>
#include <vector>
#include <boost/optional.hpp>

#include "storm-conv/settings/modules/JaniExportSettings.h"
#include "storm/storage/jani/ModelFeatures.h"

namespace storm {
    namespace converter {

        struct JaniConversionOptions {
            
            JaniConversionOptions();
            JaniConversionOptions(storm::settings::modules::JaniExportSettings const& settings);
            
            /// (Automaton,Variable)-pairs that will be transformed to location variables of the respective automaton.
            std::vector<std::pair<std::string, std::string>> locationVariables;
            
            /// If set, the model will be made standard compliant (e.g. no state rewards for discrete time models)
            bool standardCompliant;
            
            /// If set, the model is transformed into a single automaton
            bool flatten;
            
            /// If given, the model will get this name
            boost::optional<std::string> modelName;
            
            /// Only these model features are allowed in the output
            storm::jani::ModelFeatures allowedModelFeatures;
            
        };
    }
}


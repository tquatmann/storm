#include "ModelExportFormat.h"

#include <filesystem>

#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/macros.h"

namespace storm {
namespace io {

ModelExportFormat getModelExportFormatFromString(std::string const& input) {
    if (input == "dot") {
        return ModelExportFormat::Dot;
    } else if (input == "drdd") {
        return ModelExportFormat::Drdd;
    } else if (input == "drn") {
        return ModelExportFormat::Drn;
    } else if (input == "json") {
        return ModelExportFormat::Json;
    } else if (input == "dmb") {
        return ModelExportFormat::Dmb;
    }
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The model export format '" << input << "' does not match any known format.");
}

std::string toString(ModelExportFormat const& input) {
    switch (input) {
        case ModelExportFormat::Dot:
            return "dot";
        case ModelExportFormat::Drdd:
            return "drdd";
        case ModelExportFormat::Drn:
            return "drn";
        case ModelExportFormat::Json:
            return "json";
        case ModelExportFormat::Dmb:
            return "dmb";
    }
    STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "Unhandled model export format.");
}

ModelExportFormat getModelExportFormatFromFileExtension(std::string const& filename) {
    std::filesystem::path path(filename);
    if (path.has_extension()) {
        return getModelExportFormatFromString(path.extension().string().substr(1));
    } else {
        return ModelExportFormat::Dmb;
    }
}

}  // namespace io
}  // namespace storm

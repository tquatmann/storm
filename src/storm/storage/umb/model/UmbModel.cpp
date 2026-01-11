#include "storm/storage/umb/model/UmbModel.h"

#include <boost/algorithm/string/join.hpp>
#include <sstream>

#include "storm/storage/umb/model/Validation.h"

namespace storm::umb {

std::string UmbModel::getShortModelInformation() const {
    std::stringstream info;
    if (index.modelData) {
        info << index.modelData->name.value_or("Unnamed_Model");
        if (index.modelData->version) {
            info << " v" << index.modelData->version.value();
        }
    } else {
        info << "Unnamed_Model";
    }
    return info.str();
}

std::string UmbModel::getModelInformation() const {
    std::stringstream info;
    info << "UMB Model " << getShortModelInformation();
    auto optionalPrint = [&info](std::string const& label, auto const& value) {
        if (value) {
            info << "\t" << label << ": " << value.value() << "\n";
        }
    };
    if (index.modelData) {
        auto const& md = index.modelData.value();
        if (md.authors && !md.authors->empty()) {
            info << "\tAuthors: " << boost::join(md.authors.value(), "; ") << "\n";
        }
        optionalPrint("Description", md.description);
        optionalPrint("Comment", md.comment);
        optionalPrint("DOI", md.doi);
        optionalPrint("URL", md.url);
    }
    if (index.fileData) {
        auto const& fd = index.fileData.value();
        if (fd.tool) {
            info << "\tCreated with tool: " << fd.tool.value();
            if (fd.toolVersion) {
                info << " v" << fd.toolVersion.value();
            }
            info << "\n";
        }
        if (fd.creationDate) {
            info << "\tCreation date: " << fd.creationDateAsString() << "\n";
        }
    }
    return info.str();
}

bool UmbModel::validate(std::ostream& errors) const {
    return storm::umb::validate(*this, errors);
}

void UmbModel::validateOrThrow() const {
    storm::umb::validateOrThrow(*this);
}

}  // namespace storm::umb

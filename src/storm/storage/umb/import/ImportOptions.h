#pragma once

namespace storm::umb {
struct ImportOptions {
    /*!
     * Controls building of choice labelings.
     * A choice labelling will only be built if this is set to true *and* the umb model contains appropriate information.
     */
    bool buildChoiceLabeling{true};

    /*!
     * Controls building of state valuations.
     * State valuations will only be built if this is set to true *and* the umb model contains appropriate information.
     */
    bool buildStateValuations{true};
};
}  // namespace storm::umb

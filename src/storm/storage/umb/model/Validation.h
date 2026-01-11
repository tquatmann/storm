#pragma once

#include <ostream>

#include "storm/storage/umb/model/UmbModelForward.h"

namespace storm::umb {
/*!
 * Validates the given UMB model and writes potential errors to the given output stream.
 * @return true if the UMB model is valid.
 */
bool validate(storm::umb::UmbModel const& umbModel, std::ostream& errors);

/*!
 * Validates the given UMB model. If it is invalid, an exception is thrown.
 */
void validateOrThrow(storm::umb::UmbModel const& umbModel);

}  // namespace storm::umb

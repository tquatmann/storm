#pragma once

#include <string>

namespace storm::io {

enum class CompressionMode { Default, None, Gzip, Xz };

/*!
 * @return The expression mode whose string representation matches the given input
 * @throws InvalidArgumentException if the input doesn't match any known compression mode
 */
CompressionMode getCompressionModeFromString(std::string const& input);

/*!
 * @return The string representation of the given input
 */
std::string toString(CompressionMode const& input);

}  // namespace storm::io

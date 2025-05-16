#pragma once

#include <optional>
#include <set>
#include <string>

namespace storm::umb {
struct ExportOptions {
    /*!
     * The type that is used for all kinds of values. Default picks the value type of the input model.
     * @note UMB allows to use different value types for e.g. probabilities and rewards.
     *       We don't support that for now in favour of a cleaner option interface.
     */
    enum class ValueType { Default, Rational, Double, DoubleInterval } valueType{ValueType::Default};
};
}  // namespace storm::umb

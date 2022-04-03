#pragma once

#include <limits>
#include <utility>
#include "storm/solver/OptimizationDirection.h"
#include "storm/utility/macros.h"

namespace storm::utility {

/*!
 * Stores and manages an extremal value
 */
template<storm::OptimizationDirection Dir, typename ValueType, typename Enable = void>
class Extremum {
   public:
    Extremum() = default;
    Extremum(ValueType const& value) : _value(value), _empty(false) {}
    Extremum(ValueType&& value) : _value(std::move(value)), _empty(false) {}
    Extremum(Extremum const&) = default;
    Extremum(Extremum&&) = default;
    Extremum& operator=(Extremum const&) = default;
    Extremum& operator=(Extremum&&) = default;
    ~Extremum() = default;

    Extremum& operator=(ValueType const& value) {
        _value = value;
        _empty = false;
        return *this;
    }
    Extremum& operator=(ValueType&& value) {
        _value = std::move(value);
        _empty = false;
        return *this;
    }

    bool better(ValueType const& value) const {
        if constexpr (storm::solver::minimize(Dir)) {
            return _empty || value < _value;
        } else {
            static_assert(storm::solver::maximize(Dir));
            return _empty || value > _value;
        }
    }

    bool operator&=(Extremum const& other) {
        if (other.empty()) {
            return false;
        }
        return *this &= *other;
    }

    bool operator&=(Extremum&& other) {
        if (other.empty()) {
            return false;
        }
        return *this &= std::move(*other);
    }

    bool operator&=(ValueType const& value) {
        if (better(value)) {
            _value = value;
            _empty = false;
            return true;
        }
        return false;
    }

    bool operator&=(ValueType&& value) {
        if (better(value)) {
            _value = std::move(value);
            _empty = false;
            return true;
        }
        return false;
    }

    bool empty() const {
        return _empty;
    }

    ValueType const& operator*() const {
        STORM_LOG_ASSERT(!empty(), "tried to get empty extremum.");
        return _value;
    }

    ValueType& operator*() {
        STORM_LOG_ASSERT(!empty(), "tried to get empty extremum.");
        return _value;
    }

    std::optional<ValueType> getOptionalValue() const {
        if (empty()) {
            return {};
        } else {
            return _value;
        }
    }

    void reset() {
        _empty = true;
    }

   private:
    ValueType _value;
    bool _empty{true};
};

/// Specialization if -infinity and +infinity are available
template<storm::OptimizationDirection Dir, typename ValueType>
class Extremum<Dir, ValueType, typename std::enable_if_t<std::numeric_limits<ValueType>::is_iec559>> {
   public:
    Extremum() = default;
    Extremum(ValueType const& value) : _value(value){};
    Extremum(ValueType&& value) : _value(std::move(value)){};
    Extremum(Extremum const&) = default;
    Extremum(Extremum&&) = default;
    Extremum& operator=(Extremum const&) = default;
    Extremum& operator=(Extremum&& other) = default;
    ~Extremum() = default;

    Extremum& operator=(ValueType const& value) {
        _value = value;
        return *this;
    }
    Extremum& operator=(ValueType&& value) {
        _value = std::move(value);
        return *this;
    }

    bool better(ValueType const& value) const {
        if constexpr (storm::solver::minimize(Dir)) {
            return value < _value;
        } else {
            static_assert(storm::solver::maximize(Dir));
            return value > _value;
        }
    }

    bool operator&=(Extremum const& other) {
        if (other.empty()) {
            return false;
        }
        return *this &= *other;
    }

    bool operator&=(Extremum&& other) {
        if (other.empty()) {
            return false;
        }
        return *this &= std::move(*other);
    }

    bool operator&=(ValueType const& value) {
        if (better(value)) {
            _value = value;
            return true;
        }
        return false;
    }

    bool operator&=(ValueType&& value) {
        if (better(value)) {
            _value = value;
            return true;
        }
        return false;
    }

    bool empty() const {
        return _value == baseValue();
    }

    ValueType const& operator*() const {
        STORM_LOG_ASSERT(!empty(), "tried to get empty extremum.");
        return _value;
    }

    ValueType& operator*() {
        STORM_LOG_ASSERT(!empty(), "tried to get empty extremum.");
        return _value;
    }

    std::optional<ValueType> getOptionalValue() const {
        if (empty()) {
            return {};
        } else {
            return _value;
        }
    }

    void reset() {
        _value = baseValue();
    }

   private:
    ValueType constexpr baseValue() const {
        if constexpr (storm::solver::minimize(Dir)) {
            return std::numeric_limits<ValueType>::infinity();
        } else {
            static_assert(storm::solver::maximize(Dir));
            return -std::numeric_limits<ValueType>::infinity();
        }
    }

    ValueType _value{baseValue()};
};

template<typename ValueType>
using Maximum = Extremum<storm::OptimizationDirection::Maximize, ValueType>;
template<typename ValueType>
using Minimum = Extremum<storm::OptimizationDirection::Minimize, ValueType>;

}  // namespace storm::utility
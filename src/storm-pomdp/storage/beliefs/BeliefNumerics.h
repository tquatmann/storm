#pragma once

namespace storm::pomdp::beliefs {

template<typename ValueType>
struct BeliefNumerics {
    static bool lessOrEqual(ValueType const& lhs, ValueType const& rhs);
    static bool equal(ValueType const& lhs, ValueType const& rhs);
    static bool isZero(ValueType const& val);
    static bool isOne(ValueType const& val);

    /*!
     * @return a representative value that shall be used for hashing.
     */
    static ValueType valueForHash(ValueType const& val);
};

}  // namespace storm::pomdp::beliefs
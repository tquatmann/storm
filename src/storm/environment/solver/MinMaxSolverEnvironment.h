#pragma once

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/environment/SubEnvironment.h"
#include "storm/environment/solver/SolverEnvironment.h"
#include "storm/solver/MultiplicationStyle.h"
#include "storm/solver/SolverSelectionOptions.h"

namespace storm {

class MinMaxLpSolverEnvironment;

class MinMaxSolverEnvironment {
   public:
    MinMaxSolverEnvironment();
    ~MinMaxSolverEnvironment();

    storm::solver::MinMaxMethod const& getMethod() const;
    bool const& isMethodSetFromDefault() const;
    void setMethod(storm::solver::MinMaxMethod value, bool isSetFromDefault = false);
    uint64_t const& getMaximalNumberOfIterations() const;
    void setMaximalNumberOfIterations(uint64_t value);
    storm::RationalNumber const& getPrecision() const;
    void setPrecision(storm::RationalNumber value);
    bool const& getRelativeTerminationCriterion() const;
    void setRelativeTerminationCriterion(bool value);
    storm::solver::MultiplicationStyle const& getMultiplicationStyle() const;
    void setMultiplicationStyle(storm::solver::MultiplicationStyle value);
    bool isSymmetricUpdatesSet() const;
    void setSymmetricUpdates(bool value);
    MinMaxLpSolverEnvironment const& getMinMaxLpSolverEnvironment() const;
    MinMaxLpSolverEnvironment& getMinMaxLpSolverEnvironment();

   private:
    storm::solver::MinMaxMethod minMaxMethod;
    bool methodSetFromDefault;
    uint64_t maxIterationCount;
    storm::RationalNumber precision;
    bool considerRelativeTerminationCriterion;
    storm::solver::MultiplicationStyle multiplicationStyle;
    bool symmetricUpdates;
    SubEnvironment<MinMaxLpSolverEnvironment> lpEnv;
};
}  // namespace storm

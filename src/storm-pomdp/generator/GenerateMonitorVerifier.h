#pragma once
#include "storm/models/sparse/Dtmc.h"
#include "storm/models/sparse/Pomdp.h"

namespace storm::generator {

template<typename ValueType>
class GenerateMonitorVerifier {
   public:
    struct Options {
        std::string goodLabel = "good";
        std::string acceptingLabel = "accepting";
        std::string stepPrefix = "step";
        std::string horizonLabel = "horizon";
    };
    GenerateMonitorVerifier(storm::models::sparse::Dtmc<ValueType> const& mc, storm::models::sparse::Mdp<ValueType> const& monitor, Options const& options);
    std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> createProduct();

   private:
    const storm::models::sparse::Dtmc<ValueType>& mc;
    const storm::models::sparse::Mdp<ValueType>& monitor;
    Options options;
};

}  // namespace storm::generator
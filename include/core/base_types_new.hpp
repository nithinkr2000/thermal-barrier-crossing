#ifndef MCMC_BASE_TYPES_H
#define MCMC_BASE_TYPES_H

#include <NamedType/named_type.hpp>
#include <vector>

namespace mcmc::base_types {

using Beta = fluent::NamedType<float, struct BetaTag, fluent::Arithmetic>;
using Energy = fluent::NamedType<double, struct EnergyTag, fluent::Arithmetic>;
using ReplicaIndex = fluent::NamedType<int, struct RIdxTag, fluent::Comparable>;

struct Config {
    std::vector<double> v;
};
struct Matrix {
    // At 03:51:10, 18/04/2026, this is meant to be for the covariance matrix.
    // Hence, the L matrix for Cholesky decomposition is added as a member.
    // U is for LU decomposition and it'll be L.T for Cholesky.
    std::vector<std::vector<double>> M;
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;
};

} // namespace mcmc::base_types
#endif

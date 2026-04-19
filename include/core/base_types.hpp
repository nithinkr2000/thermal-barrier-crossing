#ifndef MCMC_BASE_TYPES_HPP
#define MCMC_BASE_TYPES_HPP

#include <NamedType/named_type.hpp>
#include <vector>

namespace mcmc::base_types {

using Beta = fluent::NamedType<float, struct BetaTag, fluent::Arithmetic>;
using Energy = fluent::NamedType<double, struct EnergyTag, fluent::Arithmetic>;
using ReplicaIndex = fluent::NamedType<int, struct RIdxTag, fluent::Comparable>;

struct Config {
    std::vector<double> v;
};

template <typename T> struct Matrix {
  private:
    // At 03:51:10, 18/04/2026, this is meant to be for the covariance matrix.
    // Hence, the L matrix for Cholesky decomposition is added as a member.
    // U is for LU decomposition and it'll be L.T for Cholesky.
    std::vector<std::vector<T>> M;
    std::vector<std::vector<T>> L;
    std::vector<std::vector<T>> U;

  public:
    Matrix(std::vector<std::vector<T>> matrix) : M{matrix} {};

    std::vector<std::vector<T>> &getM() const { return M; }
    std::vector<std::vector<T>> &getL() const { return L; }
    std::vector<std::vector<T>> &getU() const { return U; }

    void setM(std::vector<std::vector<T>> M1) { M = M1; }
    void setL(std::vector<std::vector<T>> M1) { L = M1; }
    void setU(std::vector<std::vector<T>> M1) { U = M1; }
    void setM(size_t rIdx, size_t cIdx, T val) { M[rIdx][cIdx] = val; }
    void setL(size_t rIdx, size_t cIdx, T val) { L[rIdx][cIdx] = val; }
    void setU(size_t rIdx, size_t cIdx, T val) { U[rIdx][cIdx] = val; }
};

} // namespace mcmc::base_types
#endif

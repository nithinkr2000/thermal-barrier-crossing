#ifndef MCMC_HELPER_FUNCS_HPP
#define MCMC_HELPER_FUNCS_HPP

#include "../core/base_types.hpp"
#include <cassert>
#include <cstddef>
#include <math.h>

namespace mcmc::utils {

template <typename T> bool InRange(const T &x, const T &lowerlim, const T &upperlim) {
    return ((x >= lowerlim) && (x <= upperlim));
}

template <typename T> std::vector<std::vector<T>> TransposeSqMat(const std::vector<std::vector<T>> &mat) {
    auto M = mat;
    size_t n = M.size();

    for (size_t rIdx = 0; rIdx < n; ++rIdx) {
        for (size_t cIdx = rIdx + 1; cIdx < n; ++cIdx) {
            std::swap(M[rIdx][cIdx], M[cIdx][rIdx]);
        }
    }

    return M;
}

template <typename T> void CheckSymmetric(const std::vector<std::vector<T>> &M, T tolerance = static_cast<T>(1e-6)) {

    size_t dim = M.size();

    for (size_t rIdx{0}; rIdx < dim; ++rIdx) {

        assert(M[rIdx].size() == dim);

        for (size_t cIdx{rIdx + 1}; cIdx < dim; ++cIdx) {

            if (std::abs(M[rIdx][cIdx] - M[cIdx][rIdx]) >= tolerance)
                throw std::runtime_error("Not a symmetric matrix.");
        }
    }
}

template <typename T> void CholeskyDecomposition(const mcmc::base_types::Matrix<T> &matrix);
} // namespace mcmc::utils
#endif

#ifndef MCMC_HELPER_FUNCS_CPP
#define MCMC_HELPER_FUNCS_CPP

#include "../../include/utils/helper_funcs.hpp"
#include <cmath>
#include <stdexcept>

namespace mcmc::utils {

template <typename T> void CholeskyDecomposition(mcmc::base_types::Matrix<T> &matrix) {

    size_t dim{matrix.getM().size()};

    assert(dim > 0);

    CheckSymmetric(matrix.getM());

    auto &M = matrix.getM();
    matrix.setL(std::vector<std::vector<T>>(dim, std::vector<T>(dim, 0)));
    auto &L = matrix.getL();

    for (size_t rIdx{0}; rIdx < dim; ++rIdx)
        for (size_t cIdx{0}; cIdx <= rIdx; ++cIdx) {

            T sumjk_1{0};

            for (size_t idx{0}; idx < cIdx; ++idx)
                sumjk_1 += L[rIdx][idx] * L[cIdx][idx];

            T val = M[rIdx][cIdx] - sumjk_1;

            if (cIdx == rIdx) {

                if (val <= T{0})
                    throw std::runtime_error("Matrix is not Symmetric Positive Definite.");

                else
                    matrix.setL(rIdx, cIdx, std::sqrt(val));
            }

            else {
                matrix.setL(rIdx, cIdx, val / matrix.getL()[cIdx][cIdx]);
            }
        }

    matrix.setU(mcmc::utils::TransposeSqMat(matrix.getL()));
}
} // namespace mcmc::utils
#endif

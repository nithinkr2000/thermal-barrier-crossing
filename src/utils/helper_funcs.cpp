#ifndef MCMC_UTILS_CPP
#define MCMC_UTILS_CPP

#include "../../include/utils/helper_funcs.hpp"

namespace mcmc::utils {
void CholeskyDecomposition(mcmc::base_types::Matrix &matrix) {

    size_t dim{matrix.M.size()};
    matrix.L = std::vector<std::vector<double>>(dim, std::vector<double>(dim, 0.0f));

    for (size_t rIdx{0}; rIdx < dim; ++rIdx)
        for (size_t cIdx{0}; cIdx <= rIdx; ++cIdx) {
            if (cIdx == rIdx) {
                double sumjk_1{0.0f};
                for (size_t idx{0}; idx < cIdx - 1; ++idx)
                    sumjk_1 += matrix.L[rIdx][idx] * matrix.L[rIdx][idx];

                matrix.L[rIdx][cIdx] = sqrt(matrix.M[rIdx][cIdx] - sumjk_1);
            } else {
                double sumjk_1{0.0f};
                for (size_t idx{0}; idx < cIdx - 1; ++idx)
                    sumjk_1 += matrix.L[cIdx][idx] * matrix.L[idx][rIdx];

                matrix.L[rIdx][cIdx] = sqrt(matrix.M[rIdx][cIdx] - sumjk_1) / matrix.L[rIdx][rIdx];
            }
        }

    matrix.U = matrix.L;
    mcmc::utils::TransposeSqMat(matrix.M);
}
} // namespace mcmc::utils
#endif

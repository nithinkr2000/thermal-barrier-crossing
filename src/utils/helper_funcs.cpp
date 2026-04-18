#include "../../include/utils/helper_funcs.hpp"

Matrix CholeskyDecomposition(const Matrix &matrix) {

    size_t dim{matrix.size()};
    Matrix lMatrix{dim, std::vector<double>(dim, 0.0f)};

    for (size_t rIdx{0}; rIdx < dim; ++rIdx)
        for (size_t cIdx{0}; cIdx <= rIdx; ++cIdx) {
            if (cIdx == rIdx) {
                double sumjk_1{0.0f};
                for (size_t idx{0}; idx < cIdx - 1; ++idx)
                    sumjk_1 += lMatrix[rIdx][idx] * lMatrix[rIdx][idx];

                lMatrix[rIdx][cIdx] = sqrt(matrix[rIdx][cIdx] - sumjk_1);
            } else {
                double sumjk_1{0.0f};
                for (size_t idx{0}; idx < cIdx - 1; ++idx)
                    sumjk_1 += lMatrix[cIdx][idx] * lMatrix[idx][rIdx];

                lMatrix[rIdx][cIdx] = sqrt(matrix[rIdx][cIdx] - sumjk_1) / lMatrix[rIdx][rIdx];
            }
        }

    return lMatrix;
}

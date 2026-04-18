#include "kernel_params.hpp"
#include <cassert>

struct GaussianKernel : KernelParams {
    double amplitude;
    Config avg;
    Matrix cov;

    double evaluate(const Config &x) {
        double e{amplitude};

        size_t dim{avg.size()};

        // Assert that the covariance matrix and mean have the same size.
        assert(dim == cov.size());
        // Assert that the input vector and mean have the same size.
        assert(dim == x.size());

        for (size_t idx{}; idx <= dim; ++idx) {
            // TODO: invert covariance matrix
        }

        return e;
    }
};

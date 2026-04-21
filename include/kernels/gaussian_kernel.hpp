#include "../core/base_types.hpp"
#include "../utils/math_helpers.hpp"
#include "kernel_params.hpp"
#include <cassert>

struct GaussianKernel : KernelParams {
    double amplitude;
    mcmc::base_types::Config avg;
    mcmc::base_types::Matrix<double> cov;

    double evaluate(const mcmc::base_types::Config &x) {

        double e{amplitude};

        size_t dim{avg.v.size()};

        auto &M = cov.getM();
        auto &L = cov.getL();
        auto &U = cov.getU();
        auto &v = x.v;

        // Assert that the covariance matrix and mean have the same size.
        assert(dim == M.size());
        assert(dim == L.size());

        // Assert that the input vector and mean have the same size.
        assert(dim == v.size());

        mcmc::base_types::Config res(std::vector<double>(dim, 0.0f));

        for (size_t idx{0}; idx < dim; ++idx)
            res.v[idx] = mcmc::math_utils::dot(L[idx], v);

        double exponent = -exponent * mcmc::math_utils::norm(res.v) / 2;
        e *= std::exp(exponent);

        return e;
    }
};

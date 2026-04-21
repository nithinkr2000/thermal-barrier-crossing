#include "../core/matrix_operators.hpp"

namespace mcmc::math_utils {

template <mcmc::base_types::VectorLike V> auto dot(const V &v1, const V &v2) {

    // Dot product between two vector types.
    // Declare return type based on the product of the first two elements.
    using S = decltype(v1[0] * v2[0]);

    S val{0};

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    for (size_t idx{0}; idx < v1.size(); ++idx)
        val = val + (v1[idx] * v2[idx]);

    return val;
}

template <mcmc::base_types::VectorLike V> double norm(const V &v) {

    // Norm of a vector type.

    double val{0.0f};

    for (size_t idx{0}; idx < v.size(); ++idx)
        val = val + (v[idx] * v[idx]);

    return val;
}

template <mcmc::base_types::VectorLike V> V mean(const V &v1, const V &v2) {

    // Element-wise mean of two vector types.

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    V res{(v1 + v2) * 0.5};

    return res;
}

template <mcmc::base_types::VectorLike V> V exp(const V &v) {

    // Element-wise exponential of vector type.

    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::exp(v[idx]);

    return res;
}

template <mcmc::base_types::VectorLike V> V sqrt(const V &v) {

    // Element-wise square-root of vector type.

    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::sqrt(v[idx]);

    return res;
}

template <mcmc::base_types::VectorLike V> V log(const V &v) {

    // Element-wise natural log of vector type.

    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::log(v[idx]);

    return res;
}

} // namespace mcmc::math_utils

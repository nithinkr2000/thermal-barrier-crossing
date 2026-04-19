#include <stdexcept>

namespace mcmc::base_types {

template <typename V>
concept VectorLike = requires(V v, size_t i) {
    v.size();
    v[i];
};

template <VectorLike V> V operator+(const V &v1, const V &v2) {

    // Overload for element-wise addition of n-D vectors.

    if (v1.size() != v2.size())
        throw std::runtime_error("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] + v2[idx];

    return res;
}

template <VectorLike V> V operator-(const V &v1, const V &v2) {

    // Overload for element-wise multiplication of n-D vectors.

    if (v1.size() != v2.size())
        throw std::runtime_error("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] - v2[idx];

    return res;
}

template <VectorLike V> V operator*(const V &v1, const V &v2) {

    // Overload for element-wise multiplication of n-D vectors.

    if (v1.size() != v2.size())
        throw std::runtime_error("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] * v2[idx];

    return res;
}

template <VectorLike V> V operator/(const V &v1, const V &v2) {

    // Overload for element-wise multiplication of n-D vectors.

    if (v1.size() != v2.size())
        throw std::runtime_error("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] / v2[idx];

    return res;
}

template <VectorLike V, typename S>
    requires std::is_arithmetic_v<S>
V operator*(const V &v, const S &s) {

    for (size_t idx{0}; idx < v.size(); ++idx)
        v[idx] = v[idx] * s;

    return v;
}

template <VectorLike V, typename S>
    requires std::is_arithmetic_v<S>
V operator*(const S &s, const V &v) {

    return v * s;
}

} // namespace mcmc::base_types

#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace mcmc::base_types {

template <typename V>
concept VectorLike = requires(V v, size_t i) {
    // Gippidy suggested using convertible_to to prevent misuse
    // Like if I try to setup a function that calculates size() for matrices.
    // That would fail in the operator overloads unless I do this.
    { v.size() } -> std::convertible_to<size_t>;
    v[i];
};

template <VectorLike V> V operator+(const V &v1, const V &v2) {

    // Overload for element-wise addition of n-D vectors.

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] + v2[idx];

    return res;
}

template <VectorLike V> V operator-(const V &v1, const V &v2) {

    // Overload for element-wise subtraction of n-D vectors.

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] - v2[idx];

    return res;
}

template <VectorLike V> V operator*(const V &v1, const V &v2) {

    // Overload for element-wise multiplication of n-D vectors.

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] * v2[idx];

    return res;
}

template <VectorLike V> V operator/(const V &v1, const V &v2) {

    // Overload for element-wise division of n-D vectors.

    if (v1.size() != v2.size())
        throw std::invalid_argument("Size mismatch between v1 and v2.");

    V res{v1};

    for (size_t idx{0}; idx < res.size(); ++idx)
        res[idx] = res[idx] / v2[idx];

    return res;
}

template <VectorLike V> V exp(const V &v) {
    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::exp(v[idx]);

    return res;
}

template <VectorLike V> V sqrt(const V &v) {
    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::sqrt(v[idx]);

    return res;
}

template <VectorLike V> V log(const V &v) {
    V res{v};

    for (size_t idx{0}; idx < v.size(); ++idx)
        res[idx] = std::log(v[idx]);

    return res;
}

template <VectorLike V, typename S>
    requires std::is_arithmetic_v<S>
V operator*(const V &v, const S &s) {

    // Overload for element-wise scaling of vectors.

    V res{v};

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

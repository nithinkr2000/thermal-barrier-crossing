#pragma once

/**
 * @file  base.hpp
 * @brief Core named types for strong typing
 *        aliases for potential and proposal functions,
 *        trajectory struct, replica data struct.
 *
 * Named types wrap raw containers and callables to prevent silent
 * argument-order errors in function signatures. Access the underlying
 * value with .get() at call sites.
 */

#include <NamedType/named_type.hpp>
#include <cassert>
#include <limits>
#include <random>
#include <utility>
#include <vector>

/**
 * ---------------------------------------------------------------------------
 *  Custom CRTP skills: exposes std::vector interface on named container types.
 *  ---------------------------------------------------------------------------
 */

// VectorAccess allows for accessing elements of the vector
// Enables ranged for loops
template <typename T> struct VectorAccess : fluent::crtp<T, VectorAccess> {

    auto size() const { return this->underlying().get().size(); }
    bool empty() const { return this->underlying().get().empty(); }

    auto &operator[](size_t idx) { return this->underlying().get()[idx]; }
    const auto &operator[](size_t idx) const { return this->underlying().get()[idx]; }

    auto &at(size_t idx) { return this->underlying().get().at(idx); }
    const auto &at(size_t idx) const { return this->underlying().get().at(idx); }

    auto &front() { return this->underlying().get().front(); }
    const auto &front() const { return this->underlying().get().front(); }

    auto &back() { return this->underlying().get().back(); }
    const auto &back() const { return this->underlying().get().back(); }

    auto begin() { return this->underlying().get().begin(); }
    auto end() { return this->underlying().get().end(); }

    auto begin() const { return this->underlying().get().cbegin(); }
    auto end() const { return this->underlying().get().cend(); }
};

// VectorMutate allows for modification of vectors
template <typename T> struct VectorMutate : fluent::crtp<T, VectorMutate> {

    void clear() { this->underlying().get().clear(); }

    template <typename U> void push_back(U &&v) { this->underlying().get().push_back(std::forward<U>(v)); }
};

// VectorMath allows for element-wise and broad-casted arithmetic
// for 1D vectors
template <typename T> struct VectorMath : fluent::crtp<T, VectorMath> {

    // operator+ allows for element-wise addition of 1D vectors
    T operator+(const T &that) const {

        const auto &a = this->underlying().get();
        const auto &b = that.get();

        assert(a.size() == b.size());

        auto result = a;
        for (size_t idx{0}; idx <= a.size(); ++idx)
            result[idx] += b[idx];

        return T(result);
    }

    // operat- allows for element-wise subtraction of 1D vectors
    T operator-(const T &that) const {

        const auto &a = this->underliying().get();
        const auto &b = that.get();

        assert(a.size() == b.size());

        auto result = a;

        for (size_t idx{0}; idx <= b.size(); ++idx)
            result[idx] -= b[idx];

        return T(result);
    }

    // operator* allows for scalar multiplication of a double and 1D vectors
    template <typename U> T operator*(U scalar) const {
        // Assert that nothing nonsensical is passed, like a string
        static_assert(std::is_arithmetic_v<U>, "Scalar must be a numeric type");
        auto result = this->underlying().get();

        for (auto &x : result)
            x *= scalar;

        return T(result);
    }

    // operator* allows for element-wise multiplication of 1D vectors
    T operator*(const T &that) const {

        const auto &a = this->underlying().get();
        const auto &b = that.get();

        assert(a.size() == b.size());
        auto result = a;

        for (size_t idx{0}; idx <= b.size(); ++idx)
            result[idx] *= b[idx];

        return T(result);
    }
};

// Single point in d-dimensional space
using Config = fluent::NamedType<std::vector<double>, struct ConfigTag, VectorAccess, VectorMutate, VectorMath>;

// Set of N-D particle positions.
using Ensemble = fluent::NamedType<std::vector<Config>, struct EnsTag, VectorAccess, VectorMutate>;

// Ordered sequence of potential energy values, parallel to a ConfigEnsemble.
using ConfigEnergy = fluent::NamedType<std::vector<double>, struct EnergyTag, VectorAccess, VectorMutate, VectorMath>;

// Per-replica inverse temperatures (β = 1/kT), one entry per exchange step.
using Betas = fluent::NamedType<std::vector<double>, struct BetasTag, VectorAccess, VectorMutate>;

// Per-replica potential-parameter index, one entry per exchange step.
using RepIdcs = fluent::NamedType<std::vector<int>, struct RepIdcsTag, VectorAccess, VectorMutate>;

/**
 * Parameters for a potential that is a sum of kernels.
 * Outer index: kernel index.
 * Inner index: parameter index within that kernel.
 * E.g. for a sum of Gaussians: params[i] = {amplitudes, means, stddevs}.
 **/
using VParams = fluent::NamedType<std::vector<std::vector<double>>, struct VParamsTag, VectorAccess>;

/**
 * Proposal function.
 * Given the current position and step size packed in a PosVec,
 * returns the proposed next position as a scalar.
 * Is Gaussian or Uniform distributed.
 */
using Proposal =
    fluent::NamedType<std::function<Ensemble(const Ensemble &, std::default_random_engine &, const Config &)>,
                      struct ProposalTag, fluent::Callable>;

// Bounds [lower, upper] confining the particle during propagation.
using Walls = fluent::NamedType<std::vector<std::pair<double, double>>, struct WallsTag, VectorAccess>;

/**
 * @brief All mutable state belonging to a single replica.
 *
 * `betas` and `repids` grow by one entry after each exchange attempt,
 * recording which inverse temperature and parameter set this replica
 * carried at each point in the trajectory.
 */
struct ReplicaState {
    Betas betas;
    RepIdcs repIdcs;
    VParams potParams;
};

struct WalkerConfig {
    Config positions;
    ConfigEnergy energies;
};

struct ReplicaInfo {
    Config x0{};
    ReplicaState state;
    WalkerConfig config;
    Walls walls{std::vector<std::pair<double, double>>{
        {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()}}};
};

/**
 * Potential energy function (surface).
 * Maps a vector of positions and a set of kernel parameters to an
 * energy value.
 */
using Potential = fluent::NamedType<std::function<ConfigEnergy(const Ensemble &, const std::vector<ReplicaInfo> &)>,
                                    struct PotFuncTag, fluent::Callable>;

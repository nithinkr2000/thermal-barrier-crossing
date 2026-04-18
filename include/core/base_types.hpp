#pragma once

/**
 * @file  base.hpp
 * @brief Core named types for strong typing.
 *
 * Named types wrap raw containers and callables to prevent silent
 * argument-order errors in function signatures. Access the underlying
 * value with .get() at call sites.
 */

#include <NamedType/named_type.hpp>
#include <NamedType/underlying_functionalities.hpp>
#include <cassert>
// #include<limits>
// #include<random>
#include <utility>
#include <vector>

/**
 * ----------------------------------------------------------------------------
 *  Custom CRTP skills: exposes std::vector interface on named container types.
 * ----------------------------------------------------------------------------
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
        for (size_t idx{0}; idx < a.size(); ++idx)
            result[idx] += b[idx];

        return T(result);
    }

    // operat- allows for element-wise subtraction of 1D vectors
    T operator-(const T &that) const {

        const auto &a = this->underlying().get();
        const auto &b = that.get();

        assert(a.size() == b.size());

        auto result = a;

        for (size_t idx{0}; idx < b.size(); ++idx)
            result[idx] -= b[idx];

        return T(result);
    }

    // operator* allows for scalar multiplication of a double and 1D vectors
    template <typename U>
        requires std::is_arithmetic_v<U>
    T operator*(U scalar) const {
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

        for (size_t idx{0}; idx < b.size(); ++idx)
            result[idx] *= b[idx];

        return T(result);
    }
};

// Scalar types for energy and beta, but using struct because NamedTypes feel like overkill
struct Energy {
    double e;
};

struct Beta {
    double b;
};

// Single point in d-dimensional space
using Config = fluent::NamedType<std::vector<double>, struct ConfigTag, VectorAccess, VectorMutate, VectorMath>;

// Set of N-D particle positions.
using ConfigEnsemble = fluent::NamedType<std::vector<Config>, struct ConfigEnsTag, VectorAccess, VectorMutate>;

// Per-replica potential energy values, parallel to a ConfigEnsemble.
using EnergyEnsemble =
    fluent::NamedType<std::vector<Energy>, struct EnergiesTag, VectorAccess, VectorMutate, VectorMath>;

// Per-replica inverse temperatures (β = 1/kT), parallel to ConfigEnsemble.
using BetaEnsemble = fluent::NamedType<std::vector<Beta>, struct BetasTag, VectorAccess, VectorMutate>;

// Per-replica Hamiltonian index, parallel to ConfigEnsemble.
using ReplicaIndices = fluent::NamedType<std::vector<int>, struct RepIdcsTag, VectorAccess, VectorMutate>;

// Bounds [lower, upper] confining the particle during propagation.
using Walls = fluent::NamedType<std::vector<std::pair<Config, Config>>, struct WallsTag, VectorAccess>;

// For covariance matrix and other cases
using Matrix = std::vector<std::vector<double>>;

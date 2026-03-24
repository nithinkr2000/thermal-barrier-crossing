#pragma once

/**
 * @file  base.hpp
 * @brief Core named types, function aliases, and the ReplicaInfo aggregate.
 *
 * Named types wrap raw containers and callables to prevent silent
 * argument-order errors in function signatures. Access the underlying
 * value with .get() at call sites.
 */

#include <vector>
#include <utility>
#include <NamedType/named_type.hpp>
#include <numbers>
#include <random>

// ---------------------------------------------------------------------------
// Custom CRTP skill: exposes std::vector interface on named container types
// ---------------------------------------------------------------------------

template <typename T1>
struct VectorInterface : fluent::crtp<T1, VectorInterface>
{
    auto   size ()  const { return this->underlying().get().size();  }
    bool   empty()  const { return this->underlying().get().empty(); }
    void   clear()        { this->underlying().get().clear();        }

    template <typename U> 
    void     push_back    (U&& v) { this->underlying().get().push_back(std::forward<U>(v)); }
    
    template <typename U> 
    void     reserve      (U&& v) { this->underlying().get().reserve(v); }

    auto begin()       { return this->underlying().get().begin(); }
    auto end  ()       { return this->underlying().get().end();   }
    auto begin() const { return this->underlying().get().begin(); }
    auto end  () const { return this->underlying().get().end();   }

    auto back()        { return this->underlying().get().back();  }
    auto back()  const { return this->underlying().get().back();  }

    // Fixed: was calling this->underlying().at(index), missing .get()
    auto&       operator[](size_t i)       { return this->underlying().get()[i];    }
    const auto& operator[](size_t i) const { return this->underlying().get()[i];    }
    auto&       at(size_t i)               { return this->underlying().get().at(i); }
    const auto& at(size_t i)         const { return this->underlying().get().at(i); }

    T1          operator+(const T1& a)     {
        assert(this->underlying().get().size() == a.size());
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] += a[i];
        return T1(temp);
    }

    T1          operator+(const T1& a) const {
        assert(this->underlying().get().size() == a.size());
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] += a[i];
        return T1(temp);
    }

    T1          operator-(const T1& a) {
        assert(this->underlying().get().size() == a.size());
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] -= a[i];
        return T1(temp);
    }

    T1          operator-(const T1& a) const {
        assert(this->underlying().get().size() == a.size());
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] -= a[i];
        return T1(temp);
    }

    template <typename T2>
    T1          operator+(const T2& a)     {
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] += a;
        return T1(temp);
    }

    template <typename T2>
    T1          operator-(const T2& a) {
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] -= a;
        return T1(temp);
    }

    template <typename T2>
    T1          operator-(const T2& a) const {
        auto temp = this->underlying().get();  // copy of the underlying vector
        for (size_t i = 0; i < temp.size(); ++i)
            temp[i] -= a;
        return T1(temp);
    }

};


// Single point in d-dimensional space
using Position = fluent::NamedType<
    std::vector<double>,
    struct PosTag,
    VectorInterface
    >;


/// Ordered sequence of N-D particle positions.
using PosVec = fluent::NamedType<
    std::vector<Position>, 
    struct PosVecTag, 
    VectorInterface>;

/// Ordered sequence of potential energy values, parallel to a PosVec.
using EVec = fluent::NamedType<
    std::vector<double>,
    struct EVecTag,
    VectorInterface>;

/// Per-replica inverse temperatures (β = 1/kT), one entry per exchange step.
using Betas = fluent::NamedType<
    std::vector<double>,
    struct BetasTag,
    VectorInterface>;

/// Per-replica potential-parameter index, one entry per exchange step.
using RepIdcs = fluent::NamedType<
    std::vector<int>,
    struct RepIdcsTag,
    VectorInterface>;

/**
 * Parameters for a potential that is a sum of kernels.
 * Outer index: kernel index. Inner index: parameter index within that kernel.
 * E.g. for a sum of Gaussians: params[i] = {amplitudes, means, stddevs}.
 */

using MultiFuncParams = fluent::NamedType<
    std::vector<std::vector<double>>,
    struct MultiFuncParamsTag,
    VectorInterface>;

/// Bounds [lower, upper] confining the particle during propagation.
using Walls = fluent::NamedType<
    std::vector<double>,
    struct WallsTag,
    VectorInterface>;


/**
 * Proposal function.
 * Given the current position and step size packed in a PosVec,
 * returns the proposed next position as a scalar.
 * Is Gaussian or Uniform distributed.
 */
using PropFunc = fluent::NamedType<
    std::function<PosVec (const PosVec&, std::default_random_engine&, double)>,
    struct PropFuncTag,
    fluent::Callable>;


/// Walls for all dimensions of the vector space
using Walls = fluent::NamedType<
    std::vector<std::pair<double, double>>, 
    struct WallsTag,
    VectorInterface>;

/**
 * @brief All mutable state belonging to a single replica.
 *
 * `betas` and `repids` grow by one entry after each exchange attempt,
 * recording which inverse temperature and parameter set this replica
 * carried at each point in the trajectory.
 */
struct ReplicaInfo
{
    Position        x0{};             ///< Current (and initial) position.
    Betas           betas{};          ///< Inverse-temperature history.
    RepIdcs         repids{};         ///< Parameter-set index history.
    MultiFuncParams vParams{};        ///< Current potential parameters.
    PosVec          positions{};      ///< Full position trajectory.
    EVec            freeEnergy{};     ///< Energy at each visited position.
    Walls           walls{            ///< Hard reflective boundaries [lo, hi].
        std::vector<std::pair<double, double>>{{
                std::numeric_limits<double>::lowest(), 
                std::numeric_limits<double>::max()
            }}};
};

/**
 * Potential energy function (surface).
 * Maps a vector of positions and a set of kernel parameters to an
 * energy value.
 */
using PotFunc = fluent::NamedType<
    std::function<EVec (const std::vector<ReplicaInfo>&)>,
    struct PotFuncTag,
    fluent::Callable>;
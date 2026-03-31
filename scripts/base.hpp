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

// ---------------------------------------------------------------------------
// Custom CRTP skill: exposes std::vector interface on named container types
// ---------------------------------------------------------------------------

template <typename T1>
struct VectorInterface : fluent::crtp<T1, VectorInterface>
{
    auto   size ()  const { return this->underlying().get().size();   }
    bool   empty()  const { return this->underlying().get().empty();  }
    void   clear()        { this->underlying().get().clear();         }

    template <typename U> void push_back (U&& v) { this->underlying().get().push_back(std::forward<U>(v)); }
    template <typename U> void reserve   (U&& v) { this->underlying().get().reserve(v);                    }

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
};


/// Ordered sequence of 1-D particle positions.
using PosVec = fluent::NamedType<
    std::vector<double>, 
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
 * E.g. for a sum of Gaussians: params[i] = {amplitude, mean, stddev}.
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


struct ProposalInput
{
    double position{0.0};  ///< Current particle position.
    double stepSize{1.0};  ///< Standard deviation or half-width of the proposal.
};


/**
 * Potential energy function.
 * Maps a vector of positions and a set of kernel parameters to a
 * parallel vector of energies.
 */
using PotFunc = fluent::NamedType<
    std::function<EVec (const PosVec&, const MultiFuncParams&)>,
    struct PotFuncTag,
    fluent::Callable>;

/**
 * Proposal function.
 * Given the current position and step size packed in a PosVec,
 * returns the proposed next position as a scalar.
 */
using PropFunc = fluent::NamedType<
    std::function<double(PosVec&, std::default_random_engine&)>,
    struct PropFuncTag,
    fluent::Callable>;


// ---------------------------------------------------------------------------
// Replica state
// ---------------------------------------------------------------------------

/**
 * @brief All mutable state belonging to a single replica.
 *
 * `betas` and `repids` grow by one entry after each exchange attempt,
 * recording which inverse temperature and parameter set this replica
 * carried at each point in the trajectory.
 */
struct ReplicaInfo
{
    double          x0{0.0};          ///< Current (and initial) position.
    Betas           betas{};          ///< Inverse-temperature history.
    RepIdcs         repids{};         ///< Parameter-set index history.
    MultiFuncParams vParams{};       ///< Potential parameters currently held.
    PosVec          positions{};      ///< Full position trajectory.
    EVec            freeEnergy{};     ///< Energy at each visited position.
    Walls           walls{            ///< Hard reflective boundaries [lo, hi].
        std::vector<double>{
        std::numeric_limits<double>::lowest(), 
        std::numeric_limits<double>::max()}};
};
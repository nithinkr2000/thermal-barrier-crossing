#ifndef KERNEL_PARAMS_HPP
#define KERNEL_PARAMS_HPP

#include "../core/base_types.hpp"
#include <NamedType/named_type.hpp>

// Abstract class for kernel functions used to construct potentials.
struct KernelParams {
    virtual ~KernelParams() = default;
    virtual double evaluate(const Config &) const = 0;
};

// using ComponentParams = std::vector<std::vector<double>>; // Rename type to improve readability
/**
 * Parameters for a potential that is a sum of kernels.
 * Outer index: kernel index.
 * Inner index: parameter index within that kernel.
 * E.g. for a sum of Gaussians: params[i] = {amplitudes, means, stddevs}.
 **/
// using VParams = fluent::NamedType<std::vector<std::vector<double>>, struct VParamsTag, VectorAccess>;

/**
 * Proposal function.
 * Given the current position and step size packed in a PosVec,
 * returns the proposed next position as a scalar.
 * Is Gaussian or Uniform distributed.
 */
/*using Proposal =
    fluent::NamedType<std::function<Ensemble(const Ensemble &, std::default_random_engine &, const Config &)>,
                      struct ProposalTag, fluent::Callable>;
*/
/**
 * @brief All mutable state belonging to a single replica.
 *
 * `betas` and `repids` grow by one entry after each exchange attempt,
 * recording which inverse temperature and parameter set this replica
 * carried at each point in the trajectory.
 */
/*struct ReplicaState {
    Betas betas;
    RepIdcs repIdcs;
    VParams potParams;
};

struct ReplicaInfo {
    Config x0{};
    ReplicaState state;
    WalkerConfig config;
    Walls walls{std::vector<std::pair<double, double>>{
        {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()}}};
};
*/
/**
 * Potential energy function (surface).
 * Maps a vector of positions and a set of kernel parameters to an
 * energy value.
 */
/*using Potential = fluent::NamedType<std::function<ConfigEnergy(const Ensemble &, const std::vector<ReplicaInfo> &)>,
                                    struct PotFuncTag, fluent::Callable>;
*/
/*
template <typename T> struct PotentialConfig {
    // Map potential function names to parameters
    // This way, after identifying and extracting parameters based on the keyword,
    // they can be stored to later reconstruct the full potential
    // Decided to use a dictionary for this
    std::unordered_map<std::string, T> component_sets;
    std::vector<StateVec> lims;
};*/

#endif

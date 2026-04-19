#ifndef MCMC_TRAJ_TYPES_HPP
#define MCMC_TRAJ_TYPES_HPP

#include "base_types.hpp"

namespace mcmc::traj_types {

struct Trajectory {
    std::vector<mcmc::base_types::Config> positions;
    std::vector<mcmc::base_types::Energy> energies;
};

// TODO: Add type for tempering simulations that hold parameters, indices

} // namespace mcmc::traj_types
#endif

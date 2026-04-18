#ifndef MCMC_TRAJ_TYPES_H
#define MCMC_TRAJ_TYPES_H

#include "base_types.hpp"

namespace mcmc::traj_types {

struct Trajectory {
    std::vector<Config> positions;
    std::vector<Energy> energies;
};

} // namespace mcmc::traj_types
#endif // !MCMC_TRAJ_TYPES_H

#ifndef PROP_MCMC_HPP
#define PROP_MCMC_HPP

#include "../core/param_sets.hpp"
#include <tuple>
#include <string>

std::tuple<StateVec, StateVec> propagate_mcmc(double& s0, double step_size,
                                              ComponentParams V_params, double beta,
                                              long unsigned int n_steps, 
                                              std::string pf, std::string proposal);


#endif
#ifndef MCMC_ONED_HPP
#define MCMC_ONED_HPP

#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <omp.h>
#include <random>
#include <functional>

using PotentialFunc = std::function<StateVec (const StateVec&, const ComponentParams& )>;
using ProposalFunc = std::function<double (const StateVec&, std::default_random_engine&)>;

StateVec BoltzmannInversion(const StateVec& s, double& beta);
double GaussianProposal(StateVec& params, std::default_random_engine& gen);
double UniformProposal(StateVec& params, std::default_random_engine& gen);

#endif
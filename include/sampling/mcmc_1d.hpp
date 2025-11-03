#ifndef MCMC_ONED_HPP
#define MCMC_ONED_HPP

#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <omp.h>
#include <random>

StateVec BoltzmannInversion(const StateVec& s, double& beta);
double GaussianProposal(StateVec& params);
double UniformProposal(StateVec& params);

#endif
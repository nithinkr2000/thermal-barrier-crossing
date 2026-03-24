#ifndef PROP_MCMC_HPP
#define PROP_MCMC_HPP

#include "../core/base.hpp"
#include "../landscape/potentials.hpp"
#include "proposals.hpp"
#include <tuple>
#include <string>
#include <random>

void PropagateMCMC( std::vector<ReplicaInfo>& repInfo,
    std::vector<MultiFuncParams> allRepParams,
    Position& stepSize,
    long nSteps,
    std::default_random_engine& gen, 
    const PotFunc& potential, 
    const PropFunc& proposal);

    
std::vector<bool> MonteCarloAcceptance(const EVec& E1, const std::vector<ReplicaInfo>&);                                              

void ReplicaExchangeMain(std::vector<ReplicaInfo>& init_reps,
                         double step_size, 
                         long n_steps, 
                         long n_ex, 
                         const PotFunc& potential,
                         const PropFunc& proposal);


#endif
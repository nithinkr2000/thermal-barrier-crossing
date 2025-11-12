#ifndef PROP_MCMC_HPP
#define PROP_MCMC_HPP

#include "../core/param_sets.hpp"
#include <tuple>
#include <string>
#include <random>

std::tuple<StateVec, StateVec> PropagateMCMC(double& s0, 
                                             double step_size,
                                             ComponentParams V_params, 
                                             double beta,
                                             long n_steps,
                                             std::default_random_engine& gen, 
                                             PotentialFunc potential, 
                                             ProposalFunc proposal);

bool ReplicaExchangeAcceptance(StateVec E, StateVec beta, std::default_random_engine& gen);                                              

struct RepInfo
{
    /*
    * @struct Contains values for storing replica state when performing 
    * replica exchange
    */

    // Initial state
    double s0{0.0};

    // Replica index
    const size_t rep_idx{0};
    
    // Potential energy params - Updated after each successful exchange
    ComponentParams V_params{};
    
    // Stores temperature after each exchange attempt
    // Stores positions and energies in continuguous vectors
    StateVec betas{}, positions{}, energies{};

};


void ReplicaExchangeMain(std::vector<RepInfo>& init_reps,
                         double step_size, 
                         long n_steps, 
                         long n_ex, 
                         PotentialFunc potential,
                         ProposalFunc proposal);


#endif
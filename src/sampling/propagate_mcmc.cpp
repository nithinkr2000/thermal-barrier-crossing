#include "../../include/core/param_sets.hpp"
#include "../../include/sampling/propagate_mcmc.hpp"
#include "../landscape/potentials.cpp"
#include "proposals.cpp"
#include <cassert>
#include <string>
#include <algorithm>



std::vector<bool> MonteCarloAcceptance(const EVec& E1, 
    const std::vector<ReplicaInfo>& repInfo) 
{
	assert(E1.size() == repInfo.size());
	
    EVec betaEnergyDiff{};
    std::vector<bool> mcmcAcceptance(E1.size(), false);

	for (int repID{}; repID < E1.size(); ++repID)
	{
        betaEnergyDiff.push_back(repInfo[repID].betas.back() * (repInfo[repID].freeEnergy.back() - E1[repID]));
        mcmcAcceptance[repID] = std::min(0.0, betaEnergyDiff.back());
    }	

	return mcmcAcceptance;
}




void PropagateMCMC(std::vector<ReplicaInfo>& repInfo, 
    Position& stepSize, 
    long nSteps, 
    std::default_random_engine& rGen, 
    const PotFunc& potential,
    const PropFunc& proposal)
{  

    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
    PosVec x1{};

    for (auto& repi: repInfo)
        x1.push_back(repi.x0);

    while(nSteps > 0)
    {
        // Proposing next position
        x1 = proposal(x1, rGen, stepSize);
        
        // Calculate corresponding energy
        EVec ENew{potential(x1, repInfo)};
        
        std::vector<bool> mcmcAcceptance{MonteCarloAcceptance(ENew, repInfo)};
        
        for (size_t repID{}; repID < mcmcAcceptance.size(); ++repID)
        {
            bool accepted = mcmcAcceptance[repID];
            auto& rep = repInfo[repID];

            rep.positions.push_back(accepted ? x1[repID] : rep.positions.back());

            if (accepted) 
                rep.x0 = rep.positions.back();

            rep.betas.push_back(rep.betas.back());
            rep.repids.push_back(rep.repids.back());

            rep.freeEnergy.push_back(accepted ? ENew[repID] : rep.freeEnergy.back());
        }
        
        --nSteps;
    };
}





void ReplicaExchangeMain(std::vector<ReplicaInfo>& init_reps,
                         double step_size, 
                         long n_steps, 
                         long n_ex, 
                         const PotFunc& potential,
                         const PropFunc& proposal)
{
    /**
     * @brief   Function to propagate the replica exchange simulation.
     *          Perform MCMC simulations, exchanges at predetermined
     *          intervals and updates positions, energies and temperatures
     *          accordingly.
     * 
     * @param   init_reps   -  The set of replicas for which simulations are 
     *                         to be performed. Contains instances of RepInfo
     *                         which contains positions, energies, temperatures,
     *                         starting structures and potential energy parameters.
     * @param   stepsize    -  The step size for the MCMC simulation.
     * @param   n_steps     -  Number of steps performed between each exchange attempt.
     * @param   n_ex        -  The number of exchanges performed in total.
     * @param   potential   -  The potential energy function.
     * @param   proposal    -  The proposal function for the next step.
     * 
     * Modifies init_reps by reference.
     */
    // Initialize random number generator
    std::random_device r;
    std::default_random_engine gen(r());
        
    for(long i = 0; i < n_ex; ++i)
    {                

        // Propagate all replicas
        for(size_t j = 0; j < init_reps.size(); ++j)
        {   
            RepInfo& curr_rep = init_reps[j];
            
            // Perform propagation for one replica
            auto [positions, energies] = PropagateMCMC(curr_rep.s0, 
                                                       step_size, 
                                                       curr_rep.V_params, 
                                                       curr_rep.betas.back(), 
                                                       n_steps, 
                                                       gen, 
                                                       potential, 
                                                       proposal);
            

            // Update positions vector
            curr_rep.positions.insert(curr_rep.positions.end(), positions.begin(), positions.end());

            // Update energies vector
            curr_rep.energies.insert(curr_rep.energies.end(), energies.begin(), energies.end());

            // Update current position
            curr_rep.s0 = curr_rep.positions.back();
        }
        
        // Initialize exchange acceptance with true
        bool ex_jj1 = true;
        size_t k = std::numeric_limits<size_t>::max();

        for (size_t j = 0; j < init_reps.size() - 1; j += 2)
        {
            // Even iterations: exchange pairs (0,1), (2,3), (4,5), ...
            if (i % 2 == 0)
                size_t k = j+1;

            // Odd iterations: exchange pairs (1,2), (3,4), (5,6), ...
            else
                size_t k = (j == init_reps.size() - 1) ? 0 : j + 1;


            ex_jj1 = ReplicaExchangeAcceptance({init_reps[j].energies.back(), init_reps[k].energies.back()}, 
                                               {init_reps[j].betas.back(), init_reps[k].betas.back()}, 
                                               gen);                
            
            // For a successful exchange
            if (ex_jj1)
            {
                // Exchange potential parameters
                std::swap(init_reps[j].V_params, init_reps[k].V_params);
                

                // Exchange temperatures
                double temp_beta = init_reps[j].betas.back();
                init_reps[j].betas.push_back(init_reps[k].betas.back());
                init_reps[k].betas.push_back(temp_beta);
            }

            // For an unsuccessful exchange            
            else
            {
                // Repeat temperatures at the end
                init_reps[j].betas.push_back(init_reps[j].betas.back());
                init_reps[k].betas.push_back(init_reps[k].betas.back());
            }
        }
    }
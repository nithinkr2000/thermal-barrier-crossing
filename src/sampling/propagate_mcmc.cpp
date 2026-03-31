#include "../../include/core/param_sets.hpp"
#include "../../include/sampling/propagate_mcmc.hpp"
#include "../landscape/potentials.cpp"
#include "proposals.cpp"
#include <cassert>
#include <string>
#include <algorithm>



void PropagateMCMC(std::vector<ReplicaInfo>& repInfo, 
    Position& stepSize, 
    long nSteps, 
    std::default_random_engine& rGen, 
    PotFunc potential,
    PropFunc proposal)
{

    while(nSteps > 0)
    {
        // Proposing next position
        PosVec x1{proposal(x0, rGen, stepSize)};
        
        // Calculate corresponding energy
        EVec ENew{potential(repInfo)};
        
        
        --nSteps;
    };

    // Pick a number from a uniform distribution in [0, 1] for the acceptance criterion
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
	
    // Initialize return vectors.
    std::vector<double> positions{};
    std::vector<double> energies{};
    std::vector<double> temp;

    for(long i = 0; i < n_steps; ++i)
    {        
        // Propose the next step       
        temp = {x0, step_size};
        double s1 = proposal(temp, gen);
        
        // Calculate energy
        std::vector<double> Es{potential({x0, s1}, V_params)};

        // Calculate Boltzmann weights for current and next states 
        // respectively.  
	    std::vector<double> probs = BoltzmannInversion({Es[1] - Es[0]}, beta);
        double p_acc = std::min(1.0, probs[0]);
        
	    // Perform the Metropolis Monte Carlo step
        double rand_val = uni_dist(gen);	
        
        if(rand_val < p_acc)
            x0 = s1;
        
        positions.push_back(s0);    // Store current position        
        energies.push_back(Es[0]);  // Store current energy

        // Every 1000 steps, reseed the generator
        if (i % 1000 == 0)
        {
            std::random_device r;
            gen.seed(r());
        }
        
    }

    // Return tuple of positions and energies
}



std::vector<bool> MonteCarloAcceptance(const EVec& E1, 
    const EVec& E2, 
    Betas& invTemperature) 
{
	assert(E1.size() == E2.size());
	
    EVec betaEnergyDiff{E1 - E2};
    std::vector<bool> mcmcAcceptance(E1.size(), false);

	for (int i{}; i < E1.size(); ++i)
	{
        betaEnergyDiff[i] *= -invTemperature[i];
        mcmcAcceptance[i] = std::min(0.0, betaEnergyDiff[i]);
    }	

	return mcmcAcceptance;
}




void ReplicaExchangeMain(std::vector<RepInfo>& init_reps,
                         double step_size, 
                         long n_steps, 
                         long n_ex, 
                         PotentialFunc potential,
                         ProposalFunc proposal)
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
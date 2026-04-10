#include "../../include/core/param_sets.hpp"
#include "../../include/sampling/propagate_mcmc.hpp"
#include "../landscape/potentials.hpp"
#include "../../include/sampling/proposals.hpp"
#include <cassert>
#include <string>
#include <algorithm>


/**
 * 
 * 
 *  MH criterion for MCMC simulators (not tempering)
 */ 
std::vector<bool> MetropolisHastingsAcceptance(const EVec& E1, 
    const std::vector<ReplicaInfo>& repInfo, 
    std::default_random_engine& rGen) 
{
	assert(E1.size() == repInfo.size());
	std::uniform_real_distribution<double> uniDist(0, 1);

    EVec betaEnergyDiff{};
    std::vector<bool> mcmcAcceptance(E1.size(), false);

	for (int repID{}; repID < E1.size(); ++repID)
	{
        betaEnergyDiff.push_back(-repInfo[repID].betas.back() * (repInfo[repID].freeEnergy.back() - E1[repID]));
        float accProb = std::min(1.0, std::exp(betaEnergyDiff.back()));
        mcmcAcceptance[repID] = (accProb == 1.0) ? true : (uniDist(rGen) < accProb);
    }	
        
	return mcmcAcceptance;
}


/**   
 * MH criterion for tempering algorithms
 */
std::vector<bool> BoltzmannInversion(const EVec& E, 
	const Betas& invTemperature, 
    std::default_random_engine& rGen) 
{
	std::uniform_real_distribution<double> uniDist(0, 1);

    EVec betaEnergyDiff{};
    std::vector<bool> mcmcAcceptance(E1.size(), false);

	for (int repID{}; repID < E1.size(); ++repID)
	{
        betaEnergyDiff.push_back(-repInfo[repID].betas.back() * (repInfo[repID].freeEnergy.back() - E1[repID]));
        float accProb = std::min(1.0, std::exp(betaEnergyDiff.back()));
        mcmcAcceptance[repID] = (accProb == 1.0) ? true : (uniDist(rGen) < accProb);
    }	
        
	return mcmcAcceptance;
}


void PropagateMCMC(std::vector<ReplicaInfo>& repInfo, 
    const Position& stepSize, 
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
        
        std::vector<bool> mcmcAcceptance{MetropolisHastingsAcceptance(ENew, repInfo, rGen)};
        
        for (size_t repIdx{}; repIdx < mcmcAcceptance.size(); ++repIdx)
        {
            bool accepted = mcmcAcceptance[repIdx];
            
            auto& rep = repInfo[repIdx];

            rep.positions.push_back(accepted ? x1[repIdx] : rep.positions.back());

            if (accepted) 
                rep.x0 = rep.positions.back();

            rep.betas.push_back(rep.betas.back());
            rep.repids.push_back(rep.repids.back());
            
            rep.freeEnergy.push_back(accepted ? ENew[repIdx] : rep.freeEnergy.back());
        }
        
        --nSteps;
    };
}



void ParallelTempering(std::vector<ReplicaInfo>& repInfo,
    const Position& stepSize,
    const int& exIdx,
    std::default_random_engine& rGen, 
    const PotFunc& potential)
{
    bool exAcceptance = true;
    std::vector<size_t> swapIdcs;
    std::vector<size_t> repIdcs;
    
    for (size_t repIdx{exIdx % 2}; repIdx < repInfo.size() - 1; repIdx += 2)
    {
        repIdcs.push_back(repIdx);
        swapIdcs.push_back(repIdx + 1);
    }

    // Handle wrap-around (exchange first and last replica) for odd exchange indices
    if (exIdx % 2 == 1)
    {
        repIdcs.push_back(repInfo.size() - 1);
        swapIdcs.push_back(0);
    }  

    for (size_t pairIdx{}; pairIdx < repIdcs.size(); ++pairIdx)
    {
        ReplicaInfo cuRep = repInfo[repIdcs[pairIdx]];
        ReplicaInfo exRep = repInfo[swapIdcs[pairIdx]];

        // Determine beta * energy differences for pair to be exchanged
        // E1 * beta1, E2 * beta2, E1 * beta2, E2 * beta1
        
        EVec bE = BoltzmannInversion(EVec(std::vector<double>{cuRep.freeEnergy.back(), 
                                                                    exRep.freeEnergy.back(), 
                                                                    cuRep.freeEnergy.back(), 
                                                                    exRep.freeEnergy.back()}), 
                                    Betas(std::vector<double>{cuRep.betas.back(), 
                                                                    exRep.betas.back(), 
                                                                    exRep.betas.back(), 
                                                                    cuRep.betas.back()}));
        
        /**
        * To avoid errors from finite precision, instead of min(1, exp(-delBeta * delE))
        * min(0, -delBeta * delE) is determined. 
        */
        
        float logProb {std::min(0.0, -(bE[0] - bE[2]) - (bE[1] - bE[3]))};
        // Generate random number, check probability accept 
        // bool accProb = (logProb > 0.0) ? 1.0: (std::exp(accProb));

        if (true)                                       
            std::cout << "Yasss";
    }                                                            
};

void ReplicaExchangeMain(std::vector<ReplicaInfo>& init_reps,
                         double step_size, 
                         long n_steps, 
                         long n_ex, 
                         const PotFunc& potential,
                         const PropFunc& proposal)
{
    
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
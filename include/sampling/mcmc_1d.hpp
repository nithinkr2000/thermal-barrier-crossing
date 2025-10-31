#ifndef MCMC_ONED_HPP
#define MCMC_ONED_HPP

#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <omp.h>
#include <random>

StateVec BoltzmannInversion(const StateVec& s, float& beta);
float GaussianProposal(StateVec& params)
{
    std::random_device r;
    std::mt19937 e1(r());
    
    std::normal_distribution<float> normal_dist(params[0], params[1]);
                       
    return normal_dist(e1);
}


#endif
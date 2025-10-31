#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <omp.h>
#include <cmath>

StateVec BoltzmannInversion(const StateVec& E, float& beta)
{
    StateVec bw(E.size());
    
    // Parallelize the exponential calculation
    #pragma omp parallel for
    for(size_t i = 0; i < E.size(); ++i) {
        bw[i] = std::exp(-beta * E[i]);
    }
    
    return bw;
}


float GaussianProposal(StateVec& params)
{
    std::random_device r;
    std::mt19937 e1(r());
    
    std::normal_distribution<float> normal_dist(params[0], params[1]);
                       
    return normal_dist(e1);
}


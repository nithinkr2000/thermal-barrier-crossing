#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <../sampling/mcmc_1d.hpp>
#include <omp.h>
#include <cmath>

StateVec BoltzmannInversion(const StateVec& E, double& beta)
{
    /*
    * @brief    Perform boltzmann inversion to determine the Boltzmann weights of
    *           the passed states.
    * 
    * @param    E      -   The energies for which the inversion is to be performed.
    * @param    beta   -   The inverse of the temperature (1/kT).
    * 
    * @return   bw     -   The inverted Boltzmann weights.
    * 
    */

    StateVec bw(E.size());
    
    // Parallelize the exponential calculation
    #pragma omp parallel for
    for(size_t i = 0; i < E.size(); ++i)
        bw[i] = std::exp(-beta * E[i]);
    
    return bw;
}


double GaussianProposal(const StateVec& params)
{
    /*
    * @brief    Generate Gaussian proposal for the next step.
    *
    * @param    params  -   Contains the current position and the step size,
    *                           which are equivalent to the mean and the standard
    *                           deviation of the proposal distribution.
    * 
    * @return   float   -   The next proposed step.
    * 
    */

    std::random_device r;
    std::mt19937 gen(r());
    
    std::normal_distribution<double> normal_dist(params[0], params[1]);
                       
    return normal_dist(gen);
}

double UniformProposal(const StateVec& params)
{
    /*
    * @brief    Generate uniform proposal for the next step.
    *
    * @param    params  -   Contains the current position and the step size,
    *                       which generates points from a window of size 
    *                       2 x step size centered at the current position. 
    * 
    * @return   float   -   The next proposed step.
    * 
    */

    std::random_device r;
    std::mt19937 gen(r());

    double left_lim (params[1] - params[0]);
    double right_lim (params[1] + params[0]);

    std::uniform_real_distribution<double> uni_dist(left_lim, right_lim);
                       
    return uni_dist(gen);
}


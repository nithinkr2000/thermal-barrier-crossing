#include <../core/param_sets.hpp>
#include <../core/helper_funcs.hpp>
#include <../sampling/mcmc_1d.hpp>
#include <omp.h>
#include <cmath>

StateVec BoltzmannInversion(const StateVec& E, double& beta)
{
    /**
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


double GaussianProposal(const StateVec& params,
                        std::default_random_engine& gen)
{
    /**
    *
    * @brief Function to propose the next move based on a Gaussian
    *        distribution about the current position.
    *       
    * @param  params - Mean and standard deviation of the Gaussian       
    *                  distribution. These are essentially the 
    *                  current position and step size.
    * 
    * @param  gen    - Generator for the random number generation.
    * 
    * @return A position drawn from the Gaussian distribution 
    *         defined using the arguments passed.
    * 
    */
    
    std::normal_distribution<double> normal_dist(params[0], params[1]);                   
    return normal_dist(gen);
}

double UniformProposal(const StateVec& params, 
                       std::default_random_engine& gen)
{
    /**
    *
    * @brief  Function to propose the next move based on a uniform 
    *         distribution about the current position
    *
    * @param  params - Center and half-width of the uniform 
    *                  distribution. These are the current
    *                  position and the step size respectively.
    *   
    * @param  gen    - Generator for random number generation
    * @return A position drawn from the uniform distribution
    *         defined using the arguments passed.
    *
    */
    
    const double from_val = params[0] - params[1];
    const double to_val = params[0] + params[1];
    
    std::uniform_real_distribution<double> uni_dist(from_val, to_val);
    
    return uni_dist(gen);
}

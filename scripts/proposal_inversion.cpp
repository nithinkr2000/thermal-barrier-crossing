#pragma once
#include <vector>
#include <random>
#include <omp.h>
#include "base.hpp"
#include <functional>


std::vector<double> BoltzmannWeight(const std::vector<double>& E, double beta)
{
    /**
    *
    * @brief Calculates Boltzmann weight based on potential energy 
    *        passed. Does not scale with partition function. 
    *
    * @param  E    - Set of energy values.
    * @param  beta - 1/kT where k is the Boltzmann constant T is the
    *               temperature.
    *
    * @return bw   - Boltzmann weights for the energies in E.
    *
    */
    
    std::vector<double> boltzWeight(std::vector<double> (E.size()) );
    
    // Parallelize the exponential calculation
    #pragma omp parallel for
    for(size_t i = 0; i < E.size(); ++i) {
        boltzWeight[i] = std::exp(-beta * E[i]);
    }
    
    return boltzWeight;
}


double gaussian_proposal(std::vector<double>& params, 
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


double uniform_proposal(std::vector<double>& params, 
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

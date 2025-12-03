#pragma once
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "base.hpp"

std::vector<double> multi_gaussian_potential(const std::vector<double>& s,
                              const std::vector<std::vector<double>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Gaussian functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Amplitudes, means and variances of the
    *                  Gaussian peaks.
    *
    * @return sumGaussians - Sum of Gaussian functions for the given 
    *                        points.
    *
    */
    std::vector<double> sumGaussians(s.size(), 0.0);
    
    #pragma omp parallel for
    for(size_t j = 0; j < s.size(); ++j) 
        for(size_t i = 0; i < params.size(); ++i) 
        {
            double numerator = (s[j] - params[i][1]) * (s[j] - params[i][1]);
            double denominator = 2 * (params[i][2] * params[i][2]);
            
            double energy_contribution = params[i][0] * std::exp( - numerator / denominator);

            sumGaussians[j] += energy_contribution;
        }
    return sumGaussians;
}


std::vector<double> multi_harmonic_potential(const std::vector<double>& s, 
                                             const std::vector<std::vector<double>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Harmonic functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Amplitudes, vertices and offsets of the 
    * 		       Harmonic wells.
    *
    * @return mhp    - Sum of Harmonic functions for the given 
    *                  points.
    *
    */
    
    std::vector<double> mhp(s.size(), std::numeric_limits<double>::infinity());

        #pragma omp parallel for
        for(size_t i = 0; i < params.size(); ++i)
            for(size_t j = 0; j < s.size(); ++j)
            {
                double contribution = (params[i][0] * 
                                    (s[j]-params[i][1]) * (s[j]-params[i][1])) +
                        params[i][2];
                                    
                mhp[j] = std::min(contribution, mhp[j]);
            }

    return mhp;
}


std::vector<double> quartic_potential(const std::vector<double>& s, 
                                     const std::vector<std::vector<double>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Harmonic functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Roots of the quartic potential, scale and 
    *                  offset. 
    *
    * @return qp     - Quartic function constructde from params.
    *
    */
    
    std::vector<double> qp(s.size(), 0.0);

    // #pragma omp parallel for
    for(size_t i = 0; i < s.size(); ++i)
	    qp[i] = ((s[i] - params[0][0]) *
	            (s[i] - params[1][0]) * 
	            (s[i] - params[2][0]) *
	            (s[i] - params[3][0]) *
	            params[4][0]) +
	            params[5][0];

    return qp;
}

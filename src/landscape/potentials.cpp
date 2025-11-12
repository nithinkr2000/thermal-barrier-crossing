#include "../../include/landscape/potentials.hpp"
#include "../../include/core/helper_funcs.hpp"
#include "../../include/core/param_sets.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <execution>
// Predefined double well potentials

StateVec GaussianDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec gdw(s.size(), 0.0f);

	#pragma omp parallel for
	for(size_t j = 0; j < s.size(); ++j) 			// Iterate over the positions (1 or 2)
		for(size_t i = 0; i < params.size(); ++i)   // Iterate over Gaussian components
		{
			double contrib = params[i][0] * std::exp(-std::pow(s[j] - params[i][1], 2) / (2 * params[i][2]*  params[i][2]));
			gdw[j] += contrib;
		}

	return gdw;
}	


StateVec HarmonicDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec hdw(s.size(), std::numeric_limits<double>::infinity());

	
	#pragma omp parallel for
	for(size_t i = 0; i < params.size(); ++i)
		for(size_t j = 0; j < s.size(); ++j)
		{
			double contrib = params[i][0] * std::pow(s[i] - params[i][1], 2) + params[i][2];
			
			hdw[j] = std::min(contrib, hdw[j]);
		}

	return hdw;
}


StateVec QuarticDoubleWell(const StateVec& s, const StateVec& params)
{
	StateVec qdw(s.size(), 0.0f);

	#pragma omp parallel for
	for(size_t i = 0; i < s.size(); ++i)
	for (size_t j = 0; j < params.size(); ++j)
		qdw[i] =( (s[i] - params[0]) * (s[i] - params[1]) * (s[i] - params[2]) * (s[i] - params[3]) * params[4]) + params[5];

	return qdw;
}


// Components for constructing double well potential

StateVec Gaussian(const StateVec& s, const StateVec& params)
{
	StateVec gaussian(s.size(), 0.0f);
	
	std::transform(s.begin(), s.end(), gaussian.begin(), [&params](double x){
		return params[0] * std::exp( -std::pow(x - params[1], 2 ) / (2 * params[2] * params[2]) ); 
	});

	return gaussian;
}


StateVec Harmonic(const StateVec& s, const StateVec& params)
{
	StateVec harmonic(s.size(), 0.0f);
	
	std::transform(s.begin(), s.end(), harmonic.begin(), [&params](double x){
		return params[0] * std::pow(x - params[1], 2) + params[2] ; 
	});

	return harmonic;
}


StateVec Polynomial(const StateVec& s, const StateVec& params)
{
	StateVec polynomial(s.size(), 0.0f);

	#pragma omp parallel for
	for(size_t i = 0; i < s.size(); ++i)
		for(size_t j = 0; j < params.size(); ++j)
			polynomial[i] += params[j] * std::pow(s[i], j);
		
	return polynomial;
}


StateVec Exponential(const StateVec& s, const StateVec& params)
{
	/**
	 * @brief 	Function to calculate exponential contributions to the 
	 * 			energy landscape. Takes the for 
	 * 			f(x) = A exp(-c (x - x0) )
	 * 			where A is scale factor, c the coefficient and x0 the
	 * 			value where the exponential is centered (= A).			
	 * 
	 * @param 	params	- Contains A, x, x0 in that order for one exactly
	 * 					  exponential functions.
	 * 
	 * @return	f(s)	- The exponential calculated for s.
	 * 
	 */

	StateVec exponential(s.size(), 0.0);
	std::transform(std::execution::par, s.begin(), s.end(), exponential.begin(), [&params](double x){
		return params[0] * std::exp( -params[1] * (x - params[2]));});
	
	return exponential;
}	

enum class ProposalTypes
{
    GAUSSIAN,
    UNIFORM,
};

enum class PotentialLandscape
{
    GAUSSIAN_DW,
    QUARTIC_DW,
    HARMONIC_DW,
    CUSTOM,
};

class FullPotential
{
	public:
		FullPotential (const std::unordered_map<std::string, std::string>& all_params)
		{
			
		}
		
		// Function to construct full potential energy surface from components
		StateVec ConstructFullPotential(const StateVec& s, const PotentialConfig<StateVec>& param_config)
		{
			StateVec potential(s.size(), 0.0f);

			size_t idx = 0;
			for(const auto& potential_params : param_config.component_sets)
			{	
				std::string potential_name = potential_params.first;
				StateVec params = potential_params.second;

				StateVec potential_lims = param_config.lims[idx];

				// If the limits are two numbers, this potential function
				// is only added in that range of positions
				if(!potential_name.compare("gaussian"))
				{	
					StateVec pots = Gaussian(s, params);
					for(auto& si : s)
						if(InRange(si, potential_lims[0], potential_lims[1]))
							potential[idx] += pots[idx];
				}	
				
				else if(!potential_name.compare("polynomial"))
				{
					StateVec pots = Polynomial(s, params);
					for(auto& si : s)
						if(InRange(si, potential_lims[0], potential_lims[1]))
							potential[idx] += pots[idx];
				}

				else if(!potential_name.compare("harmonic"))
				{
					StateVec pots = Harmonic(s, params);
					for(auto& si : s)
						if(InRange(si, potential_lims[0], potential_lims[1]))
							potential[idx] += pots[idx];
				}
				
				else
					throw std::invalid_argument("Unknown choice of component function!");
					
				++idx;
			}
			return potential;
			}
};



#include "../../include/landscape/potentials.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <stdexcept>
#include <iostream>

// Predefined double well potentials

StateVec GaussianDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec gdw(s.size(), 0.0f);

	#pragma omp parallel for collapse(2)
	for(size_t i = 0; i < params.size(); ++i)   // Iterate over Gaussian components
		for(size_t j = 0; j < s.size(); ++j) // Iterate over the positions (1 or 2)
		{
			float contrib = params[i][0] * std::exp(-std::pow(s[j] - params[i][1], 2) / (2 * params[i][2]));
			gdw[j] += contrib;
		}

	return gdw;
}	


StateVec HarmonicDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec hdw(s.size(), std::numeric_limits<float>::infinity());

	#pragma omp parallel for collapse(2)
	for(size_t i = 0; i < params.size(); ++i)
		for(size_t j = 0; j < s.size(); ++j)
		{
			float contrib = params[i][0] * std::pow(s[i] - params[i][1], 2) + params[i][2];
			
			hdw[j] = std::min(contrib, hdw[j]);
		}

	return hdw;
}


StateVec QuarticDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec qdw(s.size(), 0.0f);

	#pragma omp parallel for
	for(size_t i = 0; i < s.size(); ++i)
		qdw[i] =( (s[i] - params[0][0]) * (s[i] - params[0][1]) * (s[i] - params[0][2]) * (s[i] - params[0][3]) * params[0][4]) + params[0][5];

	return qdw;
}


// Components for constructing double well potential

StateVec Gaussian(const StateVec& s, const ComponentParams& params)
{
	StateVec gaussian(s.size(), 0.0f);
	
	std::transform(s.begin(), s.end(), gaussian.begin(), [&params](float x){
		return params[0][0] * std::exp( -std::pow(x - params[0][1], 2 ) / (2 * params[0][2]) ); 
	});

	return gaussian;
}


StateVec Harmonic(const StateVec& s, const ComponentParams& params)
{
	StateVec harmonic(s.size(), 0.0f);
	
	std::transform(s.begin(), s.end(), harmonic.begin(), [&params](float x){
		return params[0][0] * std::pow(x - params[0][1], 2) + params[0][2] ; 
	});

	return harmonic;
}


StateVec Polynomial(const StateVec& s, const ComponentParams& params)
{
	StateVec polynomial(s.size(), 0.0f);

	for(size_t i = 0; i < s.size(); ++i)
		for(size_t j = 0; j < params[0].size(); ++j)
			polynomial[i] += params[0][j] * std::pow(s[i], j);
		
	return polynomial;
}


// Function to construct full potential energy surface from components
StateVec ConstructFullPotential(const StateVec& s, const PotentialConfig& param_config)
{
	StateVec potential(s.size(), 0.0f);

	size_t idx = 0;
	for(const auto& potential_params : param_config.component_sets)
	{	
		std::string potential_name = potential_params.first;
		ComponentParams params = potential_params.second;

		StateVec potential_lims = param_config.lims[idx];

		// If the limits are two numbers, this potential function
		// is only added in that range of positions
		if(potential_name.compare("gaussian"))
		{	
			StateVec pots = Gaussian(s, params);
			for(auto& si : s)
				if((s[idx] <= potential_lims[1]) && s[idx] >= (potential_lims[0]))
					potential[idx] += pots[idx];
		}	
		
		else if(potential_name.compare("polynomial"))
		{
			StateVec pots = Polynomial(s, params);
			for(auto& si : s)
				if((s[idx] <= potential_lims[1]) && s[idx] >= (potential_lims[0]))
					potential[idx] += pots[idx];
		}

		else if(potential_name.compare("harmonic"))
		{
			StateVec pots = Harmonic(s, params);
			for(auto& si : s)
				if((s[idx] <= potential_lims[1]) && s[idx] >= (potential_lims[0]))
					potential[idx] += pots[idx];
		}

		
		else
			throw std::invalid_argument("Unknown choice of component function!");
			
		++idx;
	}
	return potential;
}
#include "../../include/landscape/potentials.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <stdexcept>
#include <iostream>

StateVec GaussianDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec gdw(s.size(), 0.0f);

	#pragma omp parallel for collapse(2)
	for(size_t i = 0; i < params.size(); ++i)   // Iterate over Gaussian components
		for(size_t j = 0; j < s.size(); ++j) // Iterate over the positions (1 or 2)
		{
			float contrib = params[i][0] * std::exp(-std::pow(s[i] - params[i][1], 2) / (2 * params[i][2]));
			gdw += contrib;
		}

	return gdw;
}	


StateVec HarmonicDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec hdw(s.ize(), std::numeric_limits<float>::infinity());

	#pragma omp parallel for collapse(2)
	for(size_t i = 0; i < params.size(); ++i)
		for(size_t j = 0; j < s.size(); ++j)
		{
			float contrib = params[i][0] * std::pow(s[i] - params[i][1], 2) + params[i][2];
			
			hdw[j] = std::min(contribution, hdw[j]);
		}

	return hdw;
}


StateVec QuarticDoubleWell(const StateVec& s, const ComponentParams& params)
{
	StateVec qdw(s.size(), 0.0f);

	#pragma omp parallel for
	for(size_t i = 0; i < s.size(); ++i)
		qdw[i] =( (s[i] - p[0]) * (s[i] - p[1]) * (s[i] - p[2]) * (s[i] - p[3]) * p[4]) + p[5];

	return qdw;
}




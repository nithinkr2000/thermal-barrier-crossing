#ifndef PARAM_SETS_HPP
#define PARAM_SETS_HPP

#include <vector>
#include <string>
#include <unordered_map>

using StateVec = std::vector<double>; 			 // Rename type to improve readability 
using ComponentParams = std::vector<std::vector<double>>; // Rename type to improve readability

template <typename T> 
struct PotentialConfig
{
	// Map potential function names to parameters
	// This way, after identifying and extracting parameters based on the keyword, 
	// they can be stored to later reconstruct the full potential
	// Decided to use a dictionary for this
	std::unordered_map <std::string, T> component_sets;
	std::vector<StateVec> lims;
};

#endif

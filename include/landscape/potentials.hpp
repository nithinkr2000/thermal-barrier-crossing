#ifndef POTENTIAL_SURFACE_HPP
#define POTENTIAL_SURFACE_HPP

#include <vector> 
#include <string>
#include "../core/param_sets.hpp" // StateVec, ComponentParams, PotentialConfig (struct for custom function)
				  

// Single functions for custom potential energy surface creation
StateVec Gaussian(const StateVec& s, const ComponentParams& params);
StateVec Harmonic(const StateVec& s, const ComponentParams& params);
StateVec Polynomial(const StateVec& s, const ComponentParams& params);

// Function to create custom potential function from parameter list
StateVec ConstructFullPotential(const StateVec& s, const PotentialConfig& param_config);

// Common double well potentials (for 1D systems i.e. motion in 1D state space)
StateVec GaussianDoubleWell(const StateVec& s, const ComponentParams& params);
StateVec HarmonicDoubleWell(const StateVec& s, const ComponentParams& params);
StateVec QuarticDoubleWell(const StateVec& s, const ComponentParams& params);

#endif 

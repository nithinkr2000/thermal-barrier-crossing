#include "../../include/core/base.hpp"
#include "../../include/landscape/potentials.hpp"
#include "../../include/core/helper_funcs.hpp"
#include <cassert>
#include <omp.h>
#include <cmath>
#include <numeric>

/**
 * To avoid errors from finite precision, instead of min(1, exp(-beta * delE))
 * min(0, -beta * delE) is determined. Hence, 
 */

EVec BoltzmannInversion(const EVec& E1, const EVec& E2, Betas invTemperature){
	assert(E1.size() == E2.size());
	EVec betaEnergyDiff = E1 - E2;

	for (int i{}; i < E1.size(); ++i)
		betaEnergyDiff[i] *= invTemperature[i];

	return betaEnergyDiff;
}


EVec PotentialGaussianBasis(const PosVec& x, const MultiFuncParams& potentialParams){
	EVec potentialEnergy(std::vector<double>(0.0, x.size()));

	#pragma omp parallel for
	for (int i{}; i < x.size(); ++i)
    for (int j{}; j < potentialParams.size(); ++j) {
        Position diff = x[i] - potentialParams[j][1]; 
        double sq_norm = std::inner_product(diff.get().begin(), diff.get().end(), 
                                            diff.get().begin(), 0.0);

        potentialEnergy[i] += potentialParams[j][0] * exp(-sq_norm / (2.0 * pow(potentialParams[j][2], 2)));
    }
	
	return potentialEnergy;
}

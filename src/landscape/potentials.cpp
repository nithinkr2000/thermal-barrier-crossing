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

EVec BoltzmannInversion(const EVec& E, 
	Betas& invTemperature)
{
	EVec betaEnergy(E);

	for (int i{}; i < E.size(); ++i)
		betaEnergy[i] *= -invTemperature[i];

	return betaEnergy;
}


EVec PotentialGaussianBasis(const std::vector<ReplicaInfo>& repInfo)
{
	size_t nReps = repInfo.size();
	EVec potentialEnergy(std::vector<double>(nReps, 0.0));
	

	#pragma omp parallel for
	for (int i{}; i < nReps; ++i)
		for (int j{}; j < repInfo[i].vParams.size(); ++j) {
			Position diff = repInfo[i].x0 - repInfo[i].vParams[j][1]; 
			double sq_norm = std::inner_product(diff.get().begin(), diff.get().end(), 
												diff.get().begin(), 0.0);

			potentialEnergy[i] += repInfo[i].vParams[j][0] * exp(-sq_norm / (2.0 * pow(repInfo[i].vParams[j][2], 2)));
		}
	
	return potentialEnergy;
}

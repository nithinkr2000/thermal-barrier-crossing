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


EVec PotentialGaussianBasis(const PosVec& x1, const std::vector<ReplicaInfo>& repInfo)
{
	size_t nReps{repInfo.size()};
	EVec potentialEnergy(std::vector<double>(nReps, 0.0));

	#pragma omp parallel for
	for (int repID{}; repID < nReps; ++repID)
	{ 
		bool flag{false};
		for (int dimID{}; dimID < x1[repID].size(); ++dimID)
		
			if ((x1[repID][dimID] < repInfo[repID].walls[dimID].first) || (x1[repID][dimID] > repInfo[repID].walls[dimID].second))
			{		
				flag = true;
				break;
			}

		if (flag)
		{
			potentialEnergy[repID] = std::numeric_limits<double>::max();	
			continue;
		}

		for (int j{}; j < repInfo[repID].vParams.size(); ++j) {
			Position diff = x1[repID] - repInfo[repID].vParams[j][1]; 
			double sq_norm = std::inner_product(diff.get().begin(), diff.get().end(), 
												diff.get().begin(), 0.0);

			potentialEnergy[repID] += repInfo[repID].vParams[j][0] * exp(-sq_norm / (2.0 * pow(repInfo[repID].vParams[j][2], 2)));
		}
	}

	return potentialEnergy;
}

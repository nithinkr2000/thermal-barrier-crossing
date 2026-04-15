#include "../../include/landscape/potentials.hpp"
#include "../../include/core/base.hpp"
#include "../../include/core/helper_funcs.hpp"
#include <cassert>
#include <cmath>
#include <numeric>

ConfigEnergy GaussianBasis(const Ensemble &x1, const std::vector<ReplicaInfo> &repInfo) {

    size_t nReps{repInfo.size()};
    ConfigEnergy potentialEnergy(std::vector<double>(nReps, 0.0));

    for (size_t repID{}; repID < nReps; ++repID) {

        bool flag{false};

        for (size_t dimID{}; dimID < x1[repID].size(); ++dimID)

            if ((x1[repID][dimID] < repInfo[repID].walls[dimID].first) ||
                (x1[repID][dimID] > repInfo[repID].walls[dimID].second)) {
                flag = true;
                break;
            }

        if (flag) {
            potentialEnergy[repID] = std::numeric_limits<double>::max();
            continue;
        }

        for (size_t j{}; j < repInfo[repID].state.potParams.size(); ++j) {

            Config diff = x1[repID][j] - repInfo[repID].state.potParams[j][1];
            double sq_norm = std::inner_product(diff.get().begin(), diff.get().end(), diff.get().begin(), 0.0);

            potentialEnergy[repID] +=
                repInfo[repID].state.potParams[j][0] *
                exp(-sq_norm / (2.0 * repInfo[repID].state.potParams[j][2] * repInfo[repID].state.potParams[j][2]));
        }
    }

    return potentialEnergy;
}

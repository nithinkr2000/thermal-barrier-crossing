#pragma once

#include "base_types.hpp"

struct Trajectory {
    std::vector<Config> positions;
    std::vector<Energy> energies;
};

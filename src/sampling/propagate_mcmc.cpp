#include "../../include/core/param_sets.hpp"
#include "../../include/landscape/potentials.hpp"
#include "../../include/sampling/propagate_mcmc.hpp"
#include "../../include/sampling/mcmc_1d.hpp"
#include <string>
#include <tuple>
#include <functional>

std::tuple<StateVec, StateVec> propagate_mcmc(double& s0, 
                                              double step_size,
                                              ComponentParams V_params, 
                                              double beta,
                                              long unsigned int n_steps, 
                                              std::function<StateVec (const StateVec&, ComponentParams&)> potential_func, 
                                              std::function<double (const StateVec&)> proposal)
{
    StateVec Es = potential_func({s0}, V_params);
    StateVec pos {s0};

    std::uniform_real_distribution<double> uni_dist(0, 1);
    std::random_device r;
    std::mt19937 gen(r());

    for (long unsigned int i = 0; i < n_steps; ++i)
    {   
        StateVec s_next {pos[i]};

        // Proposal for the next step
        s_next.push_back(proposal({pos[i]}));
        
        // Calculate energy for the current 
        // and proposed positions
        StateVec E_next = potential_func(s_next, V_params);
        StateVec probs = BoltzmannInversion(E_next, beta);
        
        double acceptance{1.0};

        if (probs[0] != 0.0)
            acceptance = std::min(1.0, probs[1] / probs[0]);
        
        double rand_val = uni_dist(gen);
        if (rand_val < acceptance)
        {
            pos.push_back(s_next[1]);
            Es.push_back(E_next[1]);
        }
          
        else
        {
            pos.push_back(s_next[0]);
            Es.push_back(E_next[0]);
        }        
    }

    return std::make_tuple(Es, pos);

}
                                
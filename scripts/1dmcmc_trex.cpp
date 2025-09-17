#include <omp.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <random>
#include <numbers>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <boost/range/irange.hpp>

std::vector<float> multi_gaussian_potential(const std::vector<float>& s, 
                                            const std::vector<std::vector<float>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Gaussian functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Amplitudes, means and variances of the
    *                  Gaussian peaks.
    *
    * @return mgp    - Sum of Gaussian functions for the given 
    *                  points.
    *
    */
    
    std::vector<float> mgp(s.size(), 0.0f);
    
    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i < params.size(); ++i) 
        for(size_t j = 0; j < s.size(); ++j) 
        {
            float contribution = params[i][0] * 
                                 std::exp(-std::pow(s[j] - params[i][1], 2) /
                                 (2 * params[i][2]));
                                 
            mgp[j] += contribution;
        }
    
    return mgp;
}


std::vector<float> multi_harmonic_potential(const std::vector<float>& s, 
                                            const std::vector<std::vector<float>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Harmonic functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Amplitudes, vertices and offsets of the 
    * 		       Harmonic wells.
    *
    * @return mhp    - Sum of Harmonic functions for the given 
    *                  points.
    *
    */
    
    std::vector<float> mhp(s.size(), 0.0f);

    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i < params.size(); ++i)
	    for(size_t j = 0; j < s.size(); ++j)
	    {
	        float contribution = params[i][0] * 
	                             std::pow(s[j]-params[i][1], 2) +
				     params[i][2];
	                             
	        mhp[j] = std::min(contribution, mhp[j]);
	    }

    return mhp;
}


std::vector<float> quartic_potential(const std::vector<float>& s, 
                                     const std::vector<std::vector<float>>& params)
{
    /**
    *
    * @brief Compute the potential energy for a set of points `s`
    *        assuming the potential energy landscape is a sum of 
    *        Harmonic functions.
    * 
    * @param  s      - Sets of 1D points.
    * @param  params - Roots of the quartic potential, scale and 
    *                  offset. 
    *
    * @return qp     - Quartic function constructde from params.
    *
    */
    
    std::vector<float> qp(s.size(), 0.0f);

    #pragma omp parallel for
    for(size_t i = 0; i < s.size(); ++i)
	    qp[i] = ((s[i] - params[0][0]) *
	            (s[i] - params[1][0]) * 
	            (s[i] - params[3][0]) *
	            (s[i] - params[3][0]) *
	            params[4][0]) +
	            params[5][0];

    return qp;
}



std::vector<float> boltzmann_weight(const std::vector<float>& E, float beta)
{
    /**
    *
    * @brief Calculates Boltzmann weight based on potential energy 
    *        passed. Does not scale with partition function. 
    *
    * @param  E    - Set of energy values.
    * @param  beta - 1/kT where k is the Boltzmann constant T is the
    *               temperature.
    *
    * @return bw   - Boltzmann weights for the energies in E.
    *
    */
    
    std::vector<float> bw(E.size());
    
    // Parallelize the exponential calculation
    #pragma omp parallel for
    for(size_t i = 0; i < E.size(); ++i) {
        bw[i] = std::exp(-beta * E[i]);
    }
    
    return bw;
}


float gaussian_proposal(std::vector<float> params, std::default_random_engine& e1)
{
    /**
    *
    * @brief Function to propose the next move based on a Gaussian
    *        distribution about the current position.
    *       
    * @param  params - Mean and standard deviation of the Gaussian       
    *                 distribution. These are essentially the 
    *                 current position and step size.
    *
    * @return A position drawn from the Gaussian distribution 
    *         defined using the arguments passed.
    * 
    */
    
    std::normal_distribution<float> normal_dist(params[0], params[1]);
                       
    return normal_dist(e1);
}


std::tuple<std::vector<float>, std::vector<float>> propagate_mcmc(float& s0, 
                                                                  const float step_size, 
                                                                  const std::vector<std::vector<float>>& V_params, 
                                                                  const float beta, 
                                                                  const int n_steps, 
                                                                  const std::string& pf, 
                                                                  std::default_random_engine& e1)
{
    /**
    *
    * @brief Generates a realization of the Markov chain using the Metropolis
    *        Monte Carlo algorithm. Returns results as a tuple of two vectors.
    *
    * @param  s0        - Starting point.
    * @param  step_size - The size of the steps to be taken i.e. standard
    *                     deviation of the Gaussian distribution used to 
    *                     propose steps.
    * @param  V_params  - Parameters for the potential energy landscape.
    * @param  beta      - 1/kT for calculating Boltzmann weights.
    * @param  n_steps   - The number of Monte Carlo steps to be taken.
    * @param  pf        - Choice of potential function for the simulation.
    *                      Must be 'gaussian' or 'harmonic'.
    *
    * @result positions - Set of positions generated by the Monte Carlo 
    *                     steps.
    * @result energies  - Set of energies for the positions visited.
    * 
    */
    
    
    
    // Initialize return vectors.
    std::vector<float> positions = {s0};
    std::vector<float> energies;

    
    for(int i = 0; i < n_steps; ++i)
    {        
        // Propose the next step 
        float s1 = gaussian_proposal({s0, step_size}, e1);

        std::vector<float> Es;
        
        // Calculate the energy based on the passed form of 
        // the potential energy surface
	    if(!pf.compare("gaussian"))
            Es = multi_gaussian_potential({s0, s1}, V_params);
        	
	    else if(!pf.compare("harmonic"))
	        Es = multi_harmonic_potential({s0, s1}, V_params);

	    else if(!pf.compare("quartic"))
	        Es = quartic_potential({s0, s1}, V_params);
		        
	    else
	        throw std::invalid_argument("Unknown potential function: " + pf);

        // Calculate Boltzmann weights for current and next states 
        // respectively.  
	    std::vector<float> probs = boltzmann_weight(Es, beta);
        
        float acceptance_probability;
        
        if(probs[0] != 0)
            acceptance_probability = std::min(1.0f, probs[1] / probs[0]);
        
        else
            acceptance_probability = 1.0f;
            
        
        // Pick a number from a uniform distribution in [0, 1]
        std::uniform_real_distribution<float> uni_dist(0.0f, 1.0f);
	
	    // Perform the Metropolis Monte Carlo step
        float rand_val = uni_dist(e1);	
        
        if(rand_val < acceptance_probability)
            s0 = s1;
        
        
        positions.push_back(s0);    // Store current position        
        energies.push_back(Es[0]);  // Store current energy
    }

    // Return tuple of positions and energies
    return std::make_tuple(positions, energies); 
}


std::vector<float> param_string_splitter(const std::string& s, char delimiter=',')
{
    /** 
    *
    * @brief Function to convert a string containing the function parameters
    *        to a vectors containing different parameters (potential, proposal
    *        Boltzmann weight, Monte Carlo steps).
    *
    * @param  s          - String to be split
    * @param  delimiter  - The delimiter for the lines with multiple parameters. 
    *
    * @return split_vals - Values obtained by splitting the strings.
    * 
    */
    
    std::vector<float> split_vals;
    
    std::stringstream ss(s);
    std::string item;
    
    while(std::getline(ss, item, delimiter))
        if(!item.empty())
            split_vals.push_back(std::stof(item));
            
    return split_vals;
}

std::vector<std::string> simple_splitter(std::string s, const std::string& delimiter) 
{
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    
    while ((pos = s.find(delimiter)) != std::string::npos) 
    {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    
    tokens.push_back(s);

    return tokens;
}


bool accept_mcmc_exchange(float E1, float E2, float beta1, float beta2, std::default_random_engine& e1)
{

    float acceptance_probability = std::min(1.0f, std::exp(-(beta1 - beta2) * (E1 - E2)));
  
    // Pick a number from a uniform distribution in [0, 1]
    std::uniform_real_distribution<float> uni_dist(0.0f, 1.0f);

    // Perform the Metropolis Monte Carlo step and return result
    return uni_dist(e1) < acceptance_probability;
}

struct rep_info
{
    float s0{0.0f};
    int rep_idx{0};
    
    // Rewrite the potential parameters every time there is an exchange    
    std::vector<std::vector<float>> V_params;
    
    // Each time an exchange is attempted, the new beta (which may be 
    // the same as in the previous step) is appended
    std::vector<float> betas;

    std::vector<float> positions, energies;
};



std::vector<rep_info> run_trex(std::vector<rep_info>& init_reps, float step_size, long n_steps, long n_ex, std::string pf, std::default_random_engine& e1)
{
    
    for(long i=0; i<n_ex; ++i)
    {
        for(size_t j=0; j<init_reps.size(); ++j)
        {   
            rep_info& curr_rep = init_reps[j];
            
            auto [trajectory, energies] = propagate_mcmc(curr_rep.s0, step_size, curr_rep.V_params, curr_rep.betas.back(), n_steps, pf, e1);
            curr_rep.positions.insert(curr_rep.positions.end(), trajectory.begin(), trajectory.end());
            curr_rep.energies.insert(curr_rep.energies.end(), energies.begin(), energies.end());
            
            curr_rep.s0 = curr_rep.positions.back();
        }
        
        bool ex_jj1 = true;
        
        // If it is an even iteration starting from 0, then do 0-1, 2-3, ...
        // If it is an odd iteration starting from 1, then do 1-2, 3-4 ...
        if(i%2 == 0)
            for(auto j: boost::irange(0, static_cast<int>(init_reps.size()), 2))
            {
                if (j + 1 >= init_reps.size())
                    throw std::runtime_error("Exchange index out of bounds");
                    
                ex_jj1 = accept_mcmc_exchange(init_reps[j].energies.back(), init_reps[j+1].energies.back(), init_reps[j].betas.back(), init_reps[j+1].betas.back(), e1);                
                
                if (ex_jj1)
                {
                    rep_info curr_rep = init_reps[j];
                    
                    init_reps[j].V_params = init_reps[j+1].V_params;
                    init_reps[j].betas.push_back(init_reps[j+1].betas.back());
                    
                    init_reps[j+1].V_params = curr_rep.V_params;
                    init_reps[j+1].betas.push_back(curr_rep.betas.back());
                }            
                else
                {
                    init_reps[j].betas.push_back(init_reps[j].betas.back());
                    init_reps[j+1].betas.push_back(init_reps[j+1].betas.back());
                }
            }
            
        else if(i%2 == 1)
            for(auto j: boost::irange(1, static_cast<int>(init_reps.size()), 2))
            {
                int next_idx = j + 1;
                
                if (j == init_reps.size()-1)
                    next_idx = 0;
                
                if (j + 1 >= init_reps.size())
                    throw std::runtime_error("Exchange index out of bounds");
                
                ex_jj1 = accept_mcmc_exchange(init_reps[j].energies.back(), init_reps[next_idx].energies.back(), init_reps[j].betas.back(), init_reps[next_idx].betas.back(), e1);
                
               
                if (ex_jj1)
                {
                    rep_info curr_rep = init_reps[j];
                    
                    init_reps[j].V_params = init_reps[next_idx].V_params;
                    init_reps[j].betas.push_back(init_reps[next_idx].betas.back());
                    
                    init_reps[next_idx].V_params = curr_rep.V_params;
                    init_reps[next_idx].betas.push_back(curr_rep.betas.back());
                }            
                else
                {
                    init_reps[j].betas.push_back(init_reps[j].betas.back());
                    init_reps[next_idx].betas.push_back(init_reps[next_idx].betas.back());
                }
            }        
    }
    
	return init_reps;
}   







int main()
{   
    
    // Initialize random number generator
    std::random_device r;
    std::default_random_engine e1(r());    
    
    std::unordered_map<std::string, std::string> global_params;
    
    std::string global_params_filename;
    std::cout << "Enter filename for other parameters (n_steps, step_size, n_ex, pf): ";
    std::cin >> global_params_filename;
    
    std::ifstream param_file(global_params_filename);

    if(!param_file.is_open())
        throw std::runtime_error("Input file '" + global_params_filename + "' not found.");
    
    std::string line;
    
    while(std::getline(param_file, line))
    {
        if(line.empty() || line[0] == '#')
            continue;
            
        // Read the line
        std::istringstream iss(line);
        std::string key, val;
        
        // Store in the unordered map i.e. dictionary
        if(std::getline(iss, key, '=') && std::getline(iss, val))                
            global_params[key] = val;
    }
    
    float step_size;
    long n_ex, n_steps;
    std::string pf = global_params["pf"];
    step_size = std::stof(global_params["step_size"]);
    n_ex = std::stol(global_params["n_ex"]);
    n_steps = std::stol(global_params["n_steps"]);
    
    std::unordered_map<std::string, std::string> all_params;
    
    std::string filenames;
    std::cout << "Entire replica input file name (prefix): ";
    std::cin >> filenames;
    
    int n_reps;
    std::cout << "Enter number of replicas: ";
    std::cin >> n_reps;
    
    std::vector<rep_info> replicas(n_reps);
    // Names of the parameters expected for the two potential functions
    std::vector<std::string> gaussian_params = {"amplitudes", "means", "variances"};
    std::vector<std::string> harmonic_params = {"scales", "vertices", "offsets"};
    std::vector<std::string> quartic_params = {"roots_scale_offset"};
        
    for(int i=1; i<n_reps+1; ++i)
    {   
        auto param_filename = filenames + std::to_string(i);
        std::ifstream param_file(param_filename);

        if(!param_file.is_open())
            throw std::runtime_error("Input file '" + param_filename + "' not found.");
        
        std::string line;
        
        while(std::getline(param_file, line))
        {
            if(line.empty() || line[0] == '#')
                continue;
                
            // Read the line
            std::istringstream iss(line);
            std::string key, val;
            
            // Store in the unordered map i.e. dictionary
            if(std::getline(iss, key, '=') && std::getline(iss, val))                
                all_params[key] = val;
        }
        
        float s0, step_size, beta, n_steps;
        
        std::vector<std::vector<float>> V_params_t;
        
        if(!pf.compare("gaussian"))
            for(auto& param : gaussian_params)
                V_params_t.push_back(param_string_splitter(all_params[param]));

        else if (!pf.compare("harmonic"))
            for(auto& param:harmonic_params)
                V_params_t.push_back(param_string_splitter(all_params[param]));
        
        else if (!pf.compare("quartic"))
            for(auto& param:quartic_params)
                V_params_t.push_back(param_string_splitter(all_params[param]));
        
        else
            throw std::invalid_argument("Unknown potential function.");
        
                
        std::vector<std::vector<float>> V_params(V_params_t[0].size(), std::vector<float>(V_params_t.size()));
        
        for(size_t i = 0; i < V_params_t.size(); ++i)
            for(size_t j = 0; j < V_params_t[0].size(); ++j)
                V_params[j][i] = V_params_t[i][j];
                
        s0 = std::stof(all_params["s0"]);
        beta = std::stof(all_params["beta"]);
        
        rep_info replica;
        replica.s0 = s0;
        replica.betas.push_back(beta);
        replica.positions.push_back(s0);
        replica.rep_idx = i;
        
        replicas.push_back(replica);
    }
    
    std::vector<rep_info> results = run_trex(replicas, step_size, n_steps, n_ex, pf, e1);

    for(size_t i=0; i<results.size(); ++i)
    {
        std::ofstream trajectory_file("trajectory" + std::to_string(i) + ".txt");
        
        for (const float& value : results[i].positions) 
            trajectory_file << value << '\n';
        
        trajectory_file.close();

        std::ofstream energy_file("energies" + std::to_string(i) + ".txt");
        
        for (const float& value : results[i].energies) 
            energy_file << value << '\n';

        energy_file.close();
    }
    
    return 0;
}




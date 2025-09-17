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


float gaussian_proposal(std::vector<float> params)
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
    std::random_device r;
    std::mt19937 e1(r());
    
    std::normal_distribution<float> normal_dist(params[0], params[1]);
                       
    return normal_dist(e1);
}


std::tuple<std::vector<float>, std::vector<float>> propagate_mcmc(float& s0, 
                                                                  float step_size, 
                                                                  std::vector<std::vector<float>> V_params, 
                                                                  float beta, 
                                                                  int n_steps, 
                                                                  std::string pf)
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
    
    // Initialize random number generator
    std::random_device r;
    std::default_random_engine e1(r());
    
    // Initialize return vectors.
    std::vector<float> positions = {s0};
    std::vector<float> energies;

    
    for(int i = 0; i < n_steps; ++i)
    {        
        // Propose the next step 
        float s1 = gaussian_proposal({s0, step_size});

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




int main()
{   
    
    std::unordered_map<std::string, std::string> all_params;
    
    std::string param_filename;
    
    std::cout << "Entire input file name: ";
    std::cin >> param_filename;
    
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
    
    // Names of the parameters expected for the two potential functions
    std::vector<std::string> gaussian_params = {"amplitudes", "means", "variances"};
    std::vector<std::string> harmonic_params = {"scales", "vertices", "offsets"};
    std::vector<std::string> quartic_params = {"roots_scale_offset"};
    
    float s0, step_size, beta, n_steps;
    std::string pf = all_params["pf"];
    
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
    step_size = std::stof(all_params["step_size"]);
    beta = std::stof(all_params["beta"]);
    n_steps = std::stoi(all_params["n_steps"]);        
   
    std::string out_suffix = all_params["out_suffix"];

    auto [trajectory, energy] = propagate_mcmc(s0, step_size, V_params, beta, n_steps, pf);

    std::ofstream trajectory_file("trajectory_" + out_suffix + ".txt");
    
    for (const float& value : trajectory) 
        trajectory_file << value << '\n';
    
    trajectory_file.close();


    std::ofstream energy_file("energies_" + out_suffix  + ".txt");
    
    for (const float& value : energy) 
        energy_file << value << '\n';

    energy_file.close();

    return 0;
}

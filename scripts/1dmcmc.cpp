/*********************************
 * Using only basic CPP.         * 
 * 1D system with a double well  *
 * potential.                    *
 * Single particle system.	 *
 *********************************/
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
#include <numeric>
#include <tuple>

std::vector<float> multi_gaussian_potential(std::vector<float> s, std::vector<std::vector<float>> params)
{
    /********************************************************************                                                         
    * Takes in s-values.                                                *     
    * Takes in parameters mean and variance for Gaussian distributions. *
    * Returns their sum - potential energy profile.                     * 
    *********************************************************************/
    std::vector<float> mgp(s.size(), 0.0f);
    
    for(size_t i=0; i < params.size(); ++i)
    {
        std::vector<float> temp(s.size());
        
        std::transform(s.begin(), s.end(), temp.begin(), 
                       [&, i](float x) {return params[i][0] * std::exp(-std::pow(x - params[i][1], 2) / (2 * params[i][2]) ); }); 
        
        std::transform(mgp.begin(), mgp.end(), temp.begin(), mgp.begin(), std::plus<float>());
    }
    
    return mgp;
}

std::vector<float> boltzmann_state_probs(std::vector<float> E, float beta)
{
    std::vector<float> bsp(E.size());
    
    std::transform(E.begin(), E.end(), bsp.begin(), [&](float x) {return std::exp(-beta * x); });
    float Z = std::accumulate(bsp.begin(), bsp.end(), 0.0f);
    
    for (auto& elem : bsp)
        elem /= Z;
    
    return bsp;
}

float gaussian_proposal(std::vector<float> params)
{
    std::random_device r;
    std::mt19937 e1(r());
    
    std::normal_distribution<float> normal_dist(params[0], params[1]);
                       
    return normal_dist(e1);
}

// Fix 1: Specify return type explicitly
std::tuple<std::vector<float>, std::vector<float>> propagate_mcmc(std::vector<float>& s, float& s0, float step_size, std::vector<std::vector<float>> V_params, float beta, int n_steps)
{
    std::random_device r;
    std::default_random_engine e1(r());
    std::vector<float> positions = {s0};
    std::vector<float> energies;  // Changed to store single energy values

    for(int i=0; i< n_steps; ++i)
    {        
        float s1 = gaussian_proposal({s0, step_size * step_size});
        std::vector<float> Es = multi_gaussian_potential({s0, s1}, V_params);
        std::vector<float> probs = boltzmann_state_probs(Es, beta);
        
        float acceptance_probability;
        if (probs[0] != 0)
            acceptance_probability = std::min(1.0f, probs[1] / probs[0]);
        else
            acceptance_probability = 1.0f;
            
        std::uniform_real_distribution<float> uni_dist(0.0f, 1.0f);
	
        float rand_val = uni_dist(e1);	
        if (rand_val < acceptance_probability)
            s0 = s1;
        
        positions.push_back(s0);        
        energies.push_back(Es[0]);  // Store current energy
    }

    return std::make_tuple(positions, energies); 
}

int main()
{   
    std::vector<std::vector<float>> V_params;
    float step_size;
    
    std::string v_params_filename, proposal_params_filename;
    std::string line;
    
    std::cout << "Enter potential parameters filename: ";
    std::cin >> v_params_filename;
    
    std::cout << "Enter proposal parameters filename: ";
    std::cin >> proposal_params_filename;
    
    // Read potential parameters (3 values per line)
    std::ifstream v_file(v_params_filename);
    if (!v_file.is_open()) 
        throw std::runtime_error("Could not open potential parameters file");
    
    while(std::getline(v_file, line))
    {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<float> row(3);
        
        if (!(iss >> row[0] >> row[1] >> row[2])) {
            throw std::invalid_argument("Invalid format in potential parameters file");
        }
        
        V_params.push_back(row);
    }
    v_file.close();
    

    std::ifstream p_file(proposal_params_filename);
    if (!p_file.is_open())
        throw std::runtime_error("Could not open proposal parameters file");

    if (!(p_file >> step_size)) {
        throw std::invalid_argument("Could not read step size from file");
    }
    p_file.close();
    
    // Example usage (you'll need to add your actual simulation parameters)
    std::vector<float> s; // Initialize as needed
    float s0 = -0.5f; // Starting position
    std::cout << "Enter starting position: ";
    std::cin >> s0;

    float beta = 1.0f; // Temperature parameter
    std::cout << "Enter beta (1/kT): ";
    std::cin >> beta;

    int n_steps{}; // Number of MCMC steps
    std::cout << "Enter the number of iterations: ";
    std::cin >> n_steps;

    auto [trajectory, energy] = propagate_mcmc(s, s0, step_size, V_params, beta, n_steps);

    std::ofstream trajectory_file("trajectory.txt");
    for (const float& value : trajectory) 
        trajectory_file << value << '\n';
    
    trajectory_file.close();


    std::ofstream energy_file("energies.txt");
    for (const float& value : energy) {
        energy_file << value << '\n';
    }
    energy_file.close();

    return 0;
}

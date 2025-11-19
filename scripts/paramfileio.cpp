#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <functional>
#include <random>
#include "potentials.cpp"
#include "proposal_inversion.cpp"
#include "base.hpp"
#include <type_traits>



std::vector<double> param_string_splitter(const std::string& s, char delimiter=',')
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
    
    std::vector<double> split_vals;
    
    std::stringstream ss(s);
    std::string item;
    
    while(std::getline(ss, item, delimiter))
        if(!item.empty())
            split_vals.push_back(std::stod(item));
            
    return split_vals;
}






std::unordered_map<std::string, std::string> read_inputfiles(std::string filename)
{

    /**
     * @brief   Function to read the input files and store raw flags
     *          as key - value pairs.
     * 
     * @param   filename    - Name of the file to be read
     * 
     * @return  all_params  - All key-val input pairs in the input file
     */

    std::unordered_map<std::string, std::string> all_params;
    std::ifstream param_file(filename);

    if(!param_file.is_open())
        throw std::runtime_error("Input file '" + filename + "' not found.");
    
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
    return all_params;
}









std::tuple<std::vector<double>, std::unordered_map<std::string, double>, std::vector<std::vector<double>>> read_potentials(std::string filename, std::string pname)
{
    std::cout << "Reading file " << filename << " with the expected potential " << pname << std::endl;
    auto all_params = read_inputfiles(filename);

    std::unordered_map<std::string, double> temp_pos;

    temp_pos["beta"] = std::stod(all_params["beta"]);      
    temp_pos["x0"] = std::stod(all_params["x0"]);

    std::vector<double> walls{param_string_splitter(all_params["walls"])};

    // Names of the parameters expected for the two potential functions
    std::vector<std::string> gaussian_params = {"amplitudes", "means", "stddevs"};
    std::vector<std::string> harmonic_params = {"scales", "vertices", "offsets"};
    std::vector<std::string> quartic_params = {"roots, scale, offset"};

    const std::unordered_map<std::string, std::vector<std::string>> paramlist{{"gaussian", gaussian_params},
                                                                               {"harmonic", harmonic_params},
                                                                               {"quartic", quartic_params}};

    std::vector<std::vector<double>> V_params_t;
    
    for (auto& pname_pvals : paramlist)
    {
        if (!pname.compare(pname_pvals.first))
        {   
            for(auto& pval: pname_pvals.second)
                V_params_t.push_back(param_string_splitter(all_params[pval]));
            break;
        }
    }

    std::vector<std::vector<double>> V_params;
        
    for(size_t i = 0; i < V_params_t.size(); ++i)
    {
        std::vector<double> temp;
        for(size_t j = 0; j < V_params_t[0].size(); ++j)
            temp.push_back(V_params_t[i][j]);

        V_params.push_back(temp);
    }    


    return std::make_tuple(walls, temp_pos, V_params);
}











std::tuple<std::unordered_map<std::string, long>, double, std::function<double (std::vector<double>&, std::default_random_engine&)>, 
std::function<std::vector<double> (const std::vector<double>&, const std::vector<std::vector<double>>&)>> 
read_mcmcparams(std::string filename)
{   
    auto all_params = read_inputfiles(filename);


    const std::unordered_map<std::string, std::function<std::vector<double> (const std::vector<double>&, 
                                                                             const std::vector<std::vector<double>>&)>> potentials
                                                                             {{"gaussian", multi_gaussian_potential}, 
                                                                              {"harmonic", multi_harmonic_potential},
                                                                              {"quartic", quartic_potential}};                                                                                    
                                                                            
    const std::unordered_map<std::string, std::function<double (std::vector<double>&, 
                                                                std::default_random_engine&)>> proposals
                                                                {{"gaussian", gaussian_proposal},            
                                                                 {"uniform", uniform_proposal}};                                                                             

    std::unordered_map<std::string, long> paramlist;
    
    paramlist["nwx"] = std::stol(all_params["nwx"]);
    paramlist["n_ex"] = std::stol(all_params["n_ex"]);
    paramlist["n_steps"] = std::stol(all_params["n_steps"]);
    double step_size = std::stod(all_params["step_size"]);

    // Potential function
    std::string pof = all_params["pof"];
    // Proposal function
    std::string prf = all_params["prf"];

    auto pot_it = potentials.find(pof);
    if (pot_it == potentials.end()) {
        throw std::runtime_error("Unknown potential: " + pof);
    }
    auto potential = pot_it->second;

    auto prop_it = proposals.find(prf);
    if (prop_it == proposals.end()) {
        throw std::runtime_error("Unknown proposal: " + prf);
    }
    auto proposal = prop_it->second;

    
    return std::make_tuple(paramlist, step_size, proposal, potential);
}



template <typename T>
void write_file(std::vector<T>& v, std::string filename, std::string format="text", char delim=',', long nwx=1)
{
    if (format.compare("bin"))
    {
        std::ofstream vwriter(filename);

        for (size_t i = 0; i < v.size(); i += nwx)
           vwriter << v[i] << delim;
        
           vwriter.close();
    }

    else
    {
        std::ofstream vwriter(filename, std::ios::binary);
        auto dataval = v.data();
        size_t element_size = sizeof(v[0]);

        for(size_t i = 0; i < v.size(); i += nwx) {
            vwriter.write(reinterpret_cast <const char*>(&dataval[i]), element_size);
        }
        vwriter.close();
    }
}
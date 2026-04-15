#include "../../include/core/io.hpp"


StateVec ParamSplitter(const std::string& s, char delimiter)
{
    /** 
    *
    * @brief    Function to convert a string containing the function parameters
    *           to a vectors containing different parameters (potential, 
    *           proposal Boltzmann weight, Monte Carlo steps).
    *
    * @param    s          - String to be split
    * @param    delimiter  - The delimiter for the lines with multiple parameters. 
    *
    * @return   split_vals - Values obtained by splitting the strings.
    * 
    */
    
    std::vector<double> split_vals;
    
    std::stringstream ss(s);
    std::string item;
    
    while(std::getline(ss, item, delimiter))
        if(!item.empty())
            split_vals.push_back(std::stof(item));
            
    return split_vals;
}


std::unordered_map<std::string, std::string> ReadParamFile(const std::string& filename)
{
    /**
     * @brief Read parameter files and extract raw parameter text values
     * 
     * @param   filename    - Name of the file to be processed
     * 
     * @result  all_params  - A list of parameters from the filename
     */
    std::unordered_map <std::string, std::string> all_params;
    
    // Open an input file stream
    std::ifstream paramfile(filename);

    if (!paramfile.is_open())
        throw std::runtime_error("Input file " + filename + " not found in directory.");
    
    // For reading lines from the file
    std::string line;

    while (std::getline(paramfile, line))
    {
        // Exclude empty lines or comments
        if(line.empty() || line[0] == '#')
            continue;
        
        // Read line into a string stream
        std::istringstream iss(line);
        std::string key, val;


        // If the key and value where both extracted, add to all parameters
        if (std::getline(iss, key, '=') && (std::getline(iss, val)))
            all_params[key] = val;
    }

    // Return
    return all_params;
}

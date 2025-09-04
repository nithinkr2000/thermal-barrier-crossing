#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>

// Currently assumed to interact via a harmonic potential of k 
// given by interaction_strength. Unique for each species-pair.


struct SpeciesProperties{
    std::string species_name;
    double species_type;
    std::map<std::string, std::vector<double>> interaction_strengths;
    // TODO: Add more properties

};

std::map<std::string, std::vector<double>> read_species_props(const std::string& filename){
    std::map<std::string, SpeciesProperties> all_species_props;
    std::ifstream file{filename};
    std::string line;
    
    while(std::getline(file, line)){
        if (line.empty() || line[0] == '#' || !isdigit(line[0])) continue;
        
        std:istringstream iss{line};
        std::string species_name;
        
        //TODO: Add a loop read the variables to the end of the line and stores their properties in a map.
        //TODO: Add checks to make sure all the species have the same number of interaction strengths
    }     
}

#ifndef IO_HPP
#define IO_HPP


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include "param_sets.hpp"

StateVec ParamSplitter(const std::string& s, char delimiter=',');
std::unordered_map<std::string, std::string> ReadParamFile(const std::string& filename);

#endif
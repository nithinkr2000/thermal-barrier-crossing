#ifndef PARAM_READER
#define PARAM_READER

#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <memory>

namespace core
{
    struct HarmonicPotential;
    struct GaussianPotential;
    struct CauchyPotential;
    
    template <typename T>
    struct is_potential : std::false_type {};
    
    template <>
    struct is_potential<HarmonicPotential> : std::true_type {};
    
    template <>
    struct is_potential<GaussianPotential> : std::true_type {};
    
    template <>
    struct is_potential<CauchyPotential> : std::true_type {};    
    
    #if __cplusplus >= 202002L
    template <typename T>
    concept Potential = is_potential<T>::value
    #endif
    
    class ParamReaderException : public std::runtime_error
    {
        public:
            // Explicit makes sure that the error message can be displayed only by calling the class explicitly
            explicit ParamReaderException(const std::string& msg) : std::runtime_error("ParamReader: " + msg) {}
            
    }    

    
    template <typename potential_type>
    class ParamReader
    {
    
        /*
            Takes filename as an argument. Checks potential type (Gaussian, Harmonic or Cauchy)
            from filename. If name does not conform to expected format, throw error saying the
            same. Otherwise, call reader. Read parameters and store them in an unordered map.
            Calculate potential energy from parameters. 
        */
        
        //  Asssertion to make sure the potential type is valid.
        static_assert(is_potential<potential_type>::value, "Parameters invalid for all implemented potentials.")
        
        public:
            
            // explicit constructor when single parameter is read - to prevent
            // implicit type conversion to other types. Needs to be string.
            explicit ParamReader(const std::string& filename);
            
            // Calculate potential from generated potential function            
            potential_type calculate_potential();
            
            // Get parameter values by name from the dictionary / unordered map
            template <typename param_type>
            param_type get_parameter(const std::string& param_name) const;
            
            
        private:
            
            // Store potential type and file name
            std::string filename_;
            std::string ptype_;
            
            // Store dictionary of parameters
            std::unordered_map <std::string, std::string> parameters_;
            
            // Internal methods to:
            // - infer potential type from name
            // - parse the passed file and get potential energy parameters
            //   which are used to calculate the potential in             
            std::string infer_potential(const std::string& filename);
            void parse_file();
            
                       
            // Using file name to detect potential type
            std::string detect_from_name(const std::string& filename) const;
            
            // Template specialization decalarations (defined after class)
            
            // Check if the parameters are as expected.
            void validate_parameters() const;
            

            // Assign potential energy parameters in `potential_type` objects. 
            // Requires `parameters_` to have been parsed and assigned by
            // private function `parse_file`. Generate a full potential function
            // from these objects that can be used for calculating 
            // energy in calculate_potential.            
            potential_type gen_potential() const;
    };
    
    struct HarmonicPotential
    {
        double vertex;
        double latus_rectum;
        
        HarmonicPotential(double v, double lr) : vertex(v), latus_rectum(lr) {}
    }
    
    struct GaussianPotential
    {
        double coef;
        double mean;
        double sigma;
        
        GaussianPotential(double c, double m, double s) : coef(c), mean(m), sigma(s) {}   
    }
    
    struct CauchyPotential
    {
        double location;
        double scale;
        
        CauchyPotential(double l, double s) : location(l), scale(s) {};
    }
    
    template <typename potential_type>
    ParamReader<potential_type>::ParamReader(const std::string& filename) : filename_(filename)
    {
        ptype_ = infer_potential(filename);
        parse_file();
        validate_parameters;
    }
    
    template <typename potential_type>
    ParamReader<potential_type>::parse_file()
    {
        std::ifstream file(filename_);
        if (!file.isopen())
            throw ParamReaderException("Cannot open " + filename_);
            
        std::string line;
        while(std::getline(file, line))
        {
            // If a line is empty or commented out
            if (line.empty() || line[0] == '#')
                continue;
            
            // Read the line
            std::istringstream iss(line);
            std::string key, val;
            
            if std::getline(iss, key, '=') && std::getline(iss, value)                
                parameters_[key] = val;
            
        }
    }
    
    template <typename potential_type>
    std::string ParamReader<potential_type>::infer_potential(const std::string& filename)
    {
        if filename.starts_with("gaussian")
            return "GaussianPotential";
        
        else if filename.starts_with("harmonic")
            return "HarmonicPotential";
            
        else if filename.starts_with("cauchy")
            return "CauchyPotential";
        
        else
            throw ParamReaderException("Cannot detect format from input filename. Change file name.");
    }
    
    
    template <typename potential_type>
    template <typename param_type>
    
    value_type ParamReader<potential_type>::get_parameter(const std::string& param_name)
    {
        auto it = parameters_.find(param_name);
        if ( it== parameters_.end())
            throw ParamReaderException("Parameter name '" + param_name + "' not found.");
            
        std::istringstream iss(it->second)
        param_type val;

        iss >> val;
        return value;
    }
    
    
    template<>
    inline void ParamReader<GaussianPotential>::validate_parameters() const 
    {
        if (ptype_ != "gaussian")
            throw ParamReaderException("Format mismatch: expected Gaussian potential energy function.")
        
        
        bool condition1 = (parameters_.find("coef") != parameters_.end() && parameters_.find("mean") != parameters_.end() && parameters_.find("sigma") != parameters_.end());
        bool condition2 = ((parameters_.find("coef").size() == parameters_.find("mean").size()) && (parameters_.find("coef").size() == parameters_.find("sigma").size()));
        
        if (!condition1)
            throw ParamReaderException("Missing parameters for Gaussian potential energy function.");

        else if (!condition2)
            throw ParamReaderException("Size mismatch in parameter sets (coefficient, mean values, standard deviations).");
            
    }
    
    template <>
    inline void ParamReader<HarmonicPotential>::validate_parameters() const
    {
        if (ptype_ != "harmonic")
            throw ParamReaderException("Format mismatch: expected Harmonic potential energy function.")
        
        
        bool condition1 = (parameters_.find("vertex") != parameters_.end() && parameters_.find("latus_rectum") != parameters_.end());
        bool condition2 = (parameters_.find("vertex").size() == parameters_.find("latus_rectum").size()) ;
        
        if (!condition1)
            throw ParamReaderException("Missing parameters for Harmonic potential energy function.");

        else if (!condition2)
            throw ParamReaderException("Size mismatch in parameter sets (vertices, latus recta).");
    }
    
    template <>
    inline void ParamReader<CauchyPotential>::validate_parameters() const
    {
        if (ptype_ != "cauchy")
            throw ParamReaderException("Format mismatch: expected Cauchy potential energy function.")
        
        
        bool condition1 = (parameters_.find("location") != parameters_.end() && parameters_.find("scale") != parameters_.end());
        bool condition2 = (parameters_.find("location").size() == parameters_.find("scale").size()) ;
        
        if (!condition1)
            throw ParamReaderException("Missing parameters for Cauchy potential energy function.");

        else if (!condition2)
            throw ParamReaderException("Size mismatch in parameter sets (location, scale).");
    }
    
    
    
    // TODO - implementing constructed potentials
    // TODO - implementing potential calculators
}


#endif









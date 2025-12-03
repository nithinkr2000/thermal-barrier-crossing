#pragma once
#include <vector>
#include <utility>
// #include <NamedType/named_type.hpp>
#include <numbers>


//template <typename T1>
//struct VectorInterface : fluent::crtp<T1, VectorInterface>
//{
//    auto size() const
//    {
//        return this->underlying().get().size();
//    }
//    
//    bool empty() const
//    {
//        return this->underlying().get().empty();
//    }
//    
//    template<typename U>
//    void push_back(U&& value)
//    {
//        this->underlying().get().push_back(std::forward<U>(value));
//    }
//
//    template <typename U_>
//    void reserve(U_&& value)
//    {
//        this->underlying().get().reserve(value);
//    }
//    
//    void clear()
//    {
//        this->underlying().get().clear();
//    }
//
//    auto begin() { return this->underlying().get().begin(); }
//    auto end() { return this->underlying().get().end(); }
//    auto begin() const { return this->underlying().get().begin(); }
//    auto end() const { return this->underlying().get().end(); }
//    auto back() { return this->underlying().get().back(); }
//    auto back() const { return this->underlying().get().back(); }
//
//    auto& at(size_t index)
//    {
//        return this->underlying().at(index);
//    }
//
//    const auto& at(size_t index) const
//    {
//        return this->underlying().at(index);
//    }
//
//    // Non-const version
//    auto& operator[](size_t index)
//    {
//        return this->underlying().get()[index];
//    }
//    
//    // Const version
//    const auto& operator[](size_t index) const
//    {
//        return this->underlying().get()[index];
//    }
//};
//
//
//// Position containers
//using PosVec = fluent::NamedType< std::vector<double>, 
//                                  struct StatePosition, 
//                                  VectorInterface,
//                                  fluent::Callable >;
//
//
//// Position containers
//using ProbVec = fluent::NamedType< std::vector<double>, 
//                                   struct StateProbability, 
//                                   VectorInterface,
//                                   fluent::Callable >;                                  
//
//
//// Energy containers
//using EVec = fluent::NamedType< std::vector<double>, 
//                                struct StateEnergy, 
//                                VectorInterface, 
//                                fluent::Callable >;
//
//
//// Inverse temperature and replica index maps
//using Betas = fluent::NamedType< std::vector< double >, 
//                                 struct InverseTempReplicaMap, 
//                                 VectorInterface,
//                                 fluent::FunctionCallable,
//                                 fluent::MethodCallable >;
//
//// Inverse temperature and replica index maps
//using RepIdcs = fluent::NamedType< std::vector< int >, 
//                                  struct InverseTempReplicaMap, 
//                                  VectorInterface,
//                                  fluent::FunctionCallable,
//                                  fluent::MethodCallable >;
//
//// Parameters for functions that are sums of multiple kernels                                      
//using MultiFuncParams = fluent::NamedType< std::vector< std::vector<double> >, 
//                                           struct SumOfFuncs, 
//                                           VectorInterface >;
//
//
//// Parameters for single functions or composing functions                                           
//using SingleFuncParams = fluent::NamedType< std::vector<double>, 
//                                            struct SingleFunc, 
//                                            VectorInterface >;
//
//
//// Potential energy functions
//using PotFunc = fluent::NamedType< std::function< EVec ( const PosVec&, 
//                                                         const MultiFuncParams& ) >, 
//                                   struct PotentialFunction, 
//                                   fluent::Callable >;
//
//
//// Proposal calculation functions                                   
//using PropFunc = fluent::NamedType< std::function< double ( PosVec&, 
//                                                            std::default_random_engine&) >, 
//                                    struct ProposalFunction, 
//                                    fluent::Callable >;
//


struct ReplicaInfo
{
    double x0{0.0};

    std::vector<int> repids{};
    std::vector<double> betas{};
    std::vector<std::vector<double>> vParams{};
    
    std::vector<double> positions{};
    std::vector<double> walls{std::numeric_limits<double>::min(), std::numeric_limits<double>::max()};
    std::vector<double> freeEnergy{};
};                                    

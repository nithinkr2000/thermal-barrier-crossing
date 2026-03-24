#include "../../include/sampling/proposals.hpp"

PosVec GaussianProposal(const PosVec& x, 
    std::default_random_engine& rGen, 
    double stepSize)
{
    PosVec nextStep(x);

    for (int repID{}; repID < nextStep.size(); ++repID)
        for (int dimID{}; dimID < nextStep[repID].size(); ++dimID){
            std::normal_distribution<double> normDist(nextStep[repID][dimID], stepSize);
            nextStep[repID][dimID] += normDist(rGen);
    }    

    return nextStep;
}


PosVec UniformProposal(const PosVec& x, 
    std::default_random_engine& rGen, 
    double stepSize)
{
    PosVec nextStep(x);

    for (int repID{}; repID < nextStep.size(); ++repID)
        for (int dimID{}; dimID < nextStep[repID].size(); ++dimID){
            
            std::uniform_real_distribution<double> uniDist(nextStep[repID][dimID] - stepSize, nextStep[repID][dimID] + stepSize);
            nextStep[repID][dimID] += uniDist(rGen);
    }    

    return nextStep;
}


std::unordered_map<std::string, PropFunc> ProposalFunction {
    {"gaussian", PropFunc(GaussianProposal)},
    {"uniform",  PropFunc(UniformProposal)}
};

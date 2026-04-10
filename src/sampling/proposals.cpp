#include "../../include/sampling/proposals.hpp"

PosVec GaussianProposal(const PosVec& x, 
    std::default_random_engine& rGen, 
    const Position& stepSize)
{
    PosVec nextStep(x);

    for (size_t repID{}; repID < nextStep.size(); ++repID)
        for (size_t dimID{}; dimID < nextStep[repID].size(); ++dimID){
            std::normal_distribution<double> normDist(nextStep[repID][dimID], stepSize[dimID]);
            nextStep[repID][dimID] = normDist(rGen);
    }    

    return nextStep;
}


PosVec UniformProposal(const PosVec& x, 
    std::default_random_engine& rGen, 
    const Position& stepSize)
{
    PosVec nextStep(x);

    for (size_t repID{}; repID < nextStep.size(); ++repID)
        for (size_t dimID{}; dimID < nextStep[repID].size(); ++dimID){
            
            std::uniform_real_distribution<double> uniDist(nextStep[repID][dimID] - stepSize[dimID], nextStep[repID][dimID] + stepSize[dimID]);
            nextStep[repID][dimID] = uniDist(rGen);
    }    

    return nextStep;
}


std::unordered_map<std::string, PropFunc> ProposalFunction {
    {"gaussian", PropFunc(GaussianProposal)},
    {"uniform",  PropFunc(UniformProposal)}
};

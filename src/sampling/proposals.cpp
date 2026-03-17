#include "../../include/sampling/proposals.hpp"

PosVec GaussianProposal(Position x, std::default_random_engine& rGen, double stepSize){
    Position nextStep(std::vector<double>(x.size(), 0.0));

    for (int i{}; i < x.size(); ++i){
        std::normal_distribution<double> normDist(x[i], stepSize);
        nextStep[i] = normDist(rGen);
    }    

    return PosVec(std::vector<Position>{x - nextStep, x + nextStep});
}


PosVec UniformProposal(Position x, std::default_random_engine& rGen, double stepSize)
{
    Position nextStep(std::vector<double>(x.size(), 0.0));

    for (int i{}; i < x.size(); ++i){
        std::uniform_real_distribution<double> uniDist(x[i] - stepSize, x[i] + stepSize);
        nextStep[i] = uniDist(rGen);
    }    

    return PosVec(std::vector<Position>{x - stepSize, x + stepSize});
}


std::unordered_map<std::string, PropFunc> ProposalFunction {
    {"gaussian", PropFunc(GaussianProposal)},
    {"uniform",  PropFunc(UniformProposal)}
};

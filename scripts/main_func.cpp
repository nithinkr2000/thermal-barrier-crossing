#include "paramfileio.cpp"
#include "1dmcmc_parallel.cpp"
#include <type_traits>
int main()
{   
    std::cout << "=== Program started ===" << std::endl;
    
    std::string param_filename;
    std::cout << "Entire input file name: ";
    std::cin >> param_filename;
    std::cout << "Got filename: " << param_filename << std::endl;

    char ex_type;
    std::cout << "Enter exchange type (h or t): ";
    std::cin >> ex_type;

    long n_reps;
    std::cout << "Enter number of replicas: " << std::flush;
    std::cin >> n_reps;
    std::cout << "Got n_reps: " << n_reps << std::endl;

    auto raw_params = read_inputfiles(param_filename + std::to_string(0));
    
    auto [all_params, step_size, proposal, potential] = read_mcmcparams(param_filename + std::to_string(0));

    long n_steps{all_params["n_steps"]};
    long n_ex{all_params["n_ex"]};
    long nwx{all_params["nwx"]};

    std::vector<ReplicaInfo> reps(n_reps);

    for(size_t i = 1; i <= n_reps; ++i)
    {
        auto [walls, beta_x0, V_params] = read_potentials(param_filename + std::to_string(i), raw_params["pof"]);
        reps[i-1].x0 = beta_x0["x0"];
        reps[i-1].betas.push_back(beta_x0["beta"]);
        reps[i-1].vParams = V_params;
        reps[i-1].walls = walls;
        
        // for (auto& params : V_params)
        // {
        //     for (auto& param: params)
        //         std::cout << param << ", ";
        //     std::cout << '\n';
        // }
        reps[i-1].repids.push_back(i);
        std::cout << "Read potential for replica " << i << "!" << std::endl;
    }

    run_trex(reps, step_size, n_steps, n_ex, potential, proposal, ex_type);
    std::cout << "=== Completed MCMC TREx" << " ===" << std::endl;

    std::cout << "=== Writing trajectory files" << " ===" << std::endl;

    std::string fname{raw_params["prefix"]};
    std::string format{raw_params["format"]};

    for(size_t i = 0; i < n_reps; ++i)
    {   
        std::string trajname{fname + "traj" + std::to_string(i + 1) + "." + format};
        write_file(reps[i].positions, trajname, format, ',', nwx);

        std::string totEname{fname + "totE" + std::to_string(i + 1) + "." + format};
        write_file(reps[i].freeEnergy, totEname, format, ',', nwx);

        std::string betaname{fname + "beta" + std::to_string(i + 1) + "." + format};
        write_file(reps[i].betas, betaname, format, ',', nwx);

        std::string idxname{fname + "idcs" + std::to_string(i + 1) + "." + format};
        write_file(reps[i].repids, idxname, format, ',', nwx);
    }
    std::cout << "=== Trajectory writing complete. Exiting program ...." << " ===" << std::endl;    
    return 0;
}

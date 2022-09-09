#include <cstdlib>
#include <iostream>

#include "Generator.hpp"
#include "CustomRandom.hpp"

int main(int argc, char * argv[])
{

    std::string conf_file = "sample.conf";
    long unsigned int seed = 1;
    int nsteps = 100;
    if (argc>1) conf_file = argv[1];
    if (argc>2) seed = std::atoi(argv[2]);
    if (argc>3) nsteps = std::atoi(argv[3]);
    std::cout << "Using conf file: " << conf_file << std::endl;
    std::cout << "Using seed: " << seed << std::endl;
    std::cout << "Running for " << nsteps << " steps.\n";

    ParamDict myParams;
    myParams.read_params(conf_file);

    std::string output_dir = "/";
    int freq = 1;
    double dt = 1e-3;
    if(myParams.is_key("output_dir")) output_dir = myParams.get_value("output_dir");
    if(myParams.is_key("freq")) freq = std::stoi(myParams.get_value("freq"));
    if(myParams.is_key("dt")) dt = std::stod(myParams.get_value("dt"));

    gsl_rng *myRNG = CustomRandom::init_rng(seed);
    Generator myNoiseGen(myParams, myRNG);

    for(int t=0; t<nsteps; t++)
    {
        std::cout << t << std::endl;
        arma::field<arma::vec> xi_r = myNoiseGen.get_xi_r();
        if(t%freq==0) myNoiseGen.save_field(xi_r, output_dir, t, dt);
        myNoiseGen.step(dt);
    }

    return 0;
}
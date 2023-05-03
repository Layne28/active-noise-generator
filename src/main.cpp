#include <cstdlib>
#include <iostream>

#include "Generator.hpp"
#include "CustomRandom.hpp"

int main(int argc, char * argv[])
{
    //Set defaults
    std::string param_file = "sample.in";
    long unsigned int seed = 1;
    int nsteps = 100;

    //Read in command line arguments if supplied
    if (argc>1) param_file = argv[1];
    if (argc>2) seed = std::atoi(argv[2]);
    if (argc>3) nsteps = std::atoi(argv[3]);

    std::cout << "Using input parameter file: " << param_file << std::endl;
    std::cout << "Using seed: " << seed << std::endl;
    std::cout << "Running for " << nsteps << " steps.\n";

    //Read in parameters from file
    ParamDict myParams;
    myParams.read_params(param_file);

    //Set parameters
    std::string output_dir = "/";
    int freq = 1;
    double dt = 1e-3;
    int do_fft = 0;
    int dim = 1;
    if(myParams.is_key("output_dir")) output_dir = myParams.get_value("output_dir");
    if(myParams.is_key("freq")) freq = std::stoi(myParams.get_value("freq"));
    if(myParams.is_key("dt")) dt = std::stod(myParams.get_value("dt"));
    if(myParams.is_key("do_fft")) do_fft = std::stoi(myParams.get_value("do_fft"));
    if(myParams.is_key("dim")) dim = std::stoi(myParams.get_value("dim"));

    //Initialize random number generator
    gsl_rng *myRNG = CustomRandom::init_rng(seed);

    //Initialize active noise generator
    //Generator3D myNoiseGen(myParams, myRNG);
    Generator myNoiseGen(myParams, myRNG);

    //Open output file for writing noise trajectory
    myNoiseGen.open_h5(output_dir);

    //Propagate an active noise trajectory
    for(int t=0; t<nsteps; t++)
    {
        if(t%100==0) std::cout << t << std::endl;
        arma::field<arma::cx_vec> xi_r = myNoiseGen.get_xi_r(do_fft);
        if(t%freq==0) myNoiseGen.save_field(xi_r, output_dir, t, dt);
        myNoiseGen.step(dt);
    }

    return 0;
}
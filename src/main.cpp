#include <cstdlib>
#include <iostream>

#include "Generator.hpp"
#include "CustomRandom.hpp"

int main(int argc, char * argv[])
{

    std::string conf_file = "sample.conf";
    long unsigned int seed = 1;
    if (argc>1) conf_file = argv[1];
    if (argc>2) seed = std::atoi(argv[2]);
    std::cout << "Using conf file: " << conf_file << std::endl;
    std::cout << "Using seed: " << seed << std::endl;

    ParamDict myParams;
    myParams.read_params(conf_file);

    gsl_rng *myRNG = CustomRandom::init_rng(seed);
    Generator myNoiseGen(myParams, myRNG);

    return 0;
}
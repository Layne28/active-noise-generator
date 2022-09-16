#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <armadillo>
#include <filesystem>
#include <math.h>
#include <complex>
#include "CustomRandom.hpp"
#include "ParamDict.hpp"

namespace fs = std::filesystem;

class Generator
{
private:

public:
    ParamDict params;
    arma::field<arma::cx_vec> xi_q;
    gsl_rng *rg;

    double dx = 1.0;
    double tau = 1.0;
    double lambda = 1.0;
    double D = 1.0;

    int nx = 4;
    int ny = 4;
    int nz = 4;

    double Lx, Ly, Lz;

    /*** Methods ***/

    //constructor
    Generator(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Generator();

    //take inverse fourier transform
    arma::field<arma::cx_vec> get_xi_r();

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier(int i, int j, int k);

    //Advance time
    void step(double dt);

    //Save noise field to file
    void save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt);
};

#endif
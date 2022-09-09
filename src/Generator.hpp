#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <armadillo>
#include <math.h>
#include <complex>
#include "CustomRandom.hpp"
#include "ParamDict.hpp"
#include "Observer.hpp"
#include "Solver.hpp"

class Generator
{
private:

public:
    ParamDict params;
    arma::field<arma::cx_vec> xi_q;
    arma::field<arma::vec> xi_r;
    gsl_rng *rg;

    double t = 0;
    double dx = 1.0;
    double tau = 1.0;
    double lambda = 1.0;
    double D = 1.0;

    int nx = 10;
    int ny = 10;
    int nz = 10;

    double Lx, Ly, Lz;

    /*** Methods ***/

    //constructor
    Generator(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Generator();

    //take inverse fourier transform
    void compute_r_from_q();

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier(int i, int j, int k);
};

#endif
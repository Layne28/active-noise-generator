#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <armadillo>
#include <experimental/filesystem>
#include <math.h>
#include <fftw3.h>
#include <complex>
#include <H5Cpp.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include "CustomRandom.hpp"
#include "ParamDict.hpp"

namespace fs = std::experimental::filesystem;

class Generator {

public:
    gsl_rng *rg;
    ParamDict params;
    arma::field<arma::cx_vec> xi_q;

    std::vector<double> spacing;
    double tau = 1.0;
    double lambda = 1.0;
    double D = 1.0;
    int dim = 3;

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
    arma::field<arma::cx_vec> get_xi_r(int do_fft=1);
    arma::field<arma::cx_vec> get_xi_r_fast();

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier_3D(int i, int j, int k);
    std::complex<double> get_rnd_gauss_fourier_2D(int i, int j);
    std::complex<double> get_rnd_gauss_fourier_1D(int i);

    //Advance time
    void step(double dt);

    //Save noise field to file
    void open_h5(std::string out_dir);
    void save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt);
};


#endif

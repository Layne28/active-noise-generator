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

//There is a parent Generator class which defines the random number generator,
//parameter dictionary, and parameters. The other classes are specific to 
//1, 2, and 3 dimensions and inherit from Generator.

class Generator {

public:    
    gsl_rng *rg;
    ParamDict params;

    double dx = 1.0;
    double tau = 1.0;
    double lambda = 1.0;
    double D = 1.0;
};

class Generator3D: public Generator {

public:
    arma::field<arma::cx_vec> xi_q;

    int nx = 4;
    int ny = 4;
    int nz = 4;

    double Lx, Ly, Lz;

    /*** Methods ***/

    //constructor
    Generator3D(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Generator3D();

    //take inverse fourier transform
    arma::field<arma::cx_vec> get_xi_r(int do_fft=0);
    arma::field<arma::cx_vec> get_xi_r_fast();

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier(int i, int j, int k);

    //Advance time
    void step(double dt);

    //Save noise field to file
    void open_h5(std::string out_dir);
    void save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt);
};

class Generator2D: public Generator {

public:
    arma::field<arma::cx_vec> xi_q;

    int nx = 4;
    int ny = 4;

    double Lx, Ly;

    /*** Methods ***/

    //constructor
    Generator2D(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Generator2D();

    //take inverse fourier transform
    arma::field<arma::cx_vec> get_xi_r(int do_fft=0);
    arma::field<arma::cx_vec> get_xi_r_fast();

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier(int i, int j);

    //Advance time
    void step(double dt);

    //Save noise field to file
    void open_h5(std::string out_dir);
    void save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt);
};

class Generator1D: public Generator {

public:
    arma::cx_vec xi_q;

    int nx = 10;

    double Lx;

    /*** Methods ***/

    //constructor
    Generator1D(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Generator1D();

    //take inverse fourier transform
    arma::cx_vec get_xi_r(int do_fft=0);
    arma::cx_vec get_xi_r_fast(); //use fftw

    //Get anticorrelated gaussian random number in fourier space
    std::complex<double> get_rnd_gauss_fourier(int i);

    //Advance time
    void step(double dt);

    //Save noise field to file
    void open_h5(std::string out_dir);
    void save_field(arma::cx_vec &theField, std::string out_dir, double t, double dt);
};

#endif

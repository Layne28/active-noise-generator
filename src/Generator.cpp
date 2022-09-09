#include "Generator.hpp"

Generator::Generator(ParamDict &theParams, gsl_rng *&the_rg)
{
    if(theParams.is_key("dx")) dx = std::stod(theParams.get_value("dx"));
    if(theParams.is_key("tau")) tau = std::stod(theParams.get_value("tau"));
    if(theParams.is_key("lambda")) lambda = std::stod(theParams.get_value("lambda"));
    if(theParams.is_key("D")) D = std::stod(theParams.get_value("D"));

    if(theParams.is_key("nx")) nx = std::stoi(theParams.get_value("nx"));
    if(theParams.is_key("ny")) ny = std::stoi(theParams.get_value("ny"));
    if(theParams.is_key("nz")) nz = std::stoi(theParams.get_value("nz"));

    Lx = nx*dx;
    Ly = ny*dx;
    Lz = nz*dx;

    xi_q = arma::field<arma::cx_vec>(nx, ny, nz);
    xi_r = arma::field<arma::vec>(nx, ny, nz);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                xi_q(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                xi_r(i,j,k) = arma::vec(3, arma::fill::zeros);
            }
        }
    }

    //set RNG
    rg = the_rg;

    //Initialize q field to ensure stationarity
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
                double prefactor = sqrt(D/tau)/(1+lambda*lambda*q_sq);
                for(int mu=0; mu<3; mu++)
                {
                    xi_q(i,j,k) = prefactor*get_rnd_gauss_fourier(i,j,k);
                }
            }
        }
    }
}

Generator::~Generator() {}

void Generator::compute_r_from_q()
{

}

std::complex<double> Generator::get_rnd_gauss_fourier(int i, int j, int k)
{
    using namespace std::complex_literals;

    if(i==0 && j==0 && k==0)
    {

    }
    else
    {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}
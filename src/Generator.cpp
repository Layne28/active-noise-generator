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
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                xi_q(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
            }
        }
    }

    //set RNG
    rg = the_rg;

    //Initialize q field to ensure stationarity
    //create a cube to track which field values have been filled
    arma::cube is_filled(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                if(is_filled((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)==1)
                {                    
                    for(int mu=0; mu<3; mu++)
                    {
                        xi_q(i,j,k)(mu) = std::conj(xi_q((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)(mu));
                    }
                    is_filled(i,j,k)=1.0;
                }
                else
                {
                    double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
                    double prefactor = sqrt(D/tau)*sqrt(Lx*Ly*Lz)/(1+lambda*lambda*q_sq);
                    for(int mu=0; mu<3; mu++)
                    {
                        xi_q(i,j,k)(mu) = prefactor*get_rnd_gauss_fourier(i,j,k);
                    }
                    is_filled(i,j,k)=1.0;
                }
            }
        }
    }
    std::cout << "Initialized generator." << std::endl;
}

Generator::~Generator() {}

arma::field<arma::cx_vec> Generator::get_xi_r(int do_fft)
{

    if(do_fft==1) return get_xi_r_fast();

    using namespace std::complex_literals;

    arma::field<arma::cx_vec> xi_r(nx, ny, nz);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                xi_r(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
            }
        }
    }

    //Not-so-fast Fourier transform
    for(int mu=0; mu<3; mu++)
    {
        for(int i=0; i<nx; i++)
        {
            for(int j=0; j<ny; j++)
            {
                for(int k=0; k<nz; k++)
                {
                    for(int q1=0; q1<nx; q1++)
                    {
                        for(int q2=0; q2<ny; q2++)
                        {
                            for(int q3=0; q3<nz; q3++)
                            {
                                std::complex<double> update = xi_q(q1,q2,q3)(mu)*std::exp(2*M_PI*1i*(
                                    (1.0*i*q1)/(1.0*nx) +
                                    (1.0*j*q2)/(1.0*ny) +
                                    (1.0*k*q3)/(1.0*nz)));
                                //Check that imaginary part is zero
                                //if(update.imag()>1e-3) std::cout << "imaginary part:" << update.imag() << std::endl;
                                xi_r(i,j,k)(mu) += update.real();
                            }
                        }
                    }
                    xi_r(i,j,k)(mu) *= 1.0/(Lx*Ly*Lz);
                }
            }
        }
    }
    return xi_r;
}

arma::field<arma::cx_vec> Generator::get_xi_r_fast()
{
    using namespace std::complex_literals;

    arma::field<arma::cx_vec> xi_r(nx, ny, nz);
    arma::cx_cube xi_rx(nx, ny, nz, arma::fill::zeros);
    arma::cx_cube xi_ry(nx, ny, nz, arma::fill::zeros);
    arma::cx_cube xi_rz(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                xi_rx(i,j,k) = xi_q(i,j,k)(0);
                xi_ry(i,j,k) = xi_q(i,j,k)(1);
                xi_rz(i,j,k) = xi_q(i,j,k)(2);
                xi_r(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                //xi_r(i,j,k)(mu) = xi_q(i,j,k)(mu);
            }
        }
    }

    //Fast Fourier transform w/ FFTW

    //std::cout << "doing fft" << std::endl;

    fftw_complex *in_x, *in_y, *in_z;
    fftw_plan px, py, pz;

    //Might be faster to initialize the plan in the constructor since all these transforms have the same size
    //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    in_x = reinterpret_cast<fftw_complex*>(xi_rx.memptr());
    in_y = reinterpret_cast<fftw_complex*>(xi_ry.memptr());
    in_z = reinterpret_cast<fftw_complex*>(xi_rz.memptr());
    px = fftw_plan_dft_3d(nx, ny, nz, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);
    py = fftw_plan_dft_3d(nx, ny, nz, in_y, in_y, FFTW_BACKWARD, FFTW_ESTIMATE);
    pz = fftw_plan_dft_3d(nx, ny, nz, in_z, in_z, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(px); 
    fftw_execute(py);
    fftw_execute(pz);

    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                xi_r(i,j,k)(0) = xi_rx(i,j,k)/(Lx*Ly*Lz);
                xi_r(i,j,k)(1) = xi_ry(i,j,k)/(Lx*Ly*Lz);
                xi_r(i,j,k)(2) = xi_rz(i,j,k)/(Lx*Ly*Lz);
            }
        }
    }
    //std::cout << xi_r << std::endl;
    fftw_destroy_plan(px);
    fftw_destroy_plan(py);
    fftw_destroy_plan(pz);
    return xi_r;
}

std::complex<double> Generator::get_rnd_gauss_fourier(int i, int j, int k)
{
    using namespace std::complex_literals;

    if( (i==0 && j==0 && k==0) ||
        (i==nx/2 && j==0 && k==0) ||
        (i==0 && j==ny/2 && k==0) ||
        (i==0 && j==0 && k==nz/2) ||
        (i==nx/2 && j==ny/2 && k==0) ||
        (i==nx/2 && j==0 && k==nz/2) ||
        (i==0 && j==ny/2 && k==nz/2) ||
        (i==nx/2 && j==ny/2 && k==nz/2) )
    {
        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else
    {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

void Generator::step(double dt)
{
    arma::field<arma::cx_vec> noise_incr(nx, ny, nz);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                noise_incr(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
            }
        }
    }

    arma::cube is_filled(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                if(is_filled((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)==1)
                {
                    for(int mu=0; mu<3; mu++)
                    {
                        noise_incr(i,j,k)(mu) = std::conj(noise_incr((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)(mu));
                    }
                    is_filled(i,j,k)=1.0;
                }
                else
                {
                    double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
                    double prefactor = sqrt(2*D*dt)*sqrt(Lx*Ly*Lz)/tau/(1+lambda*lambda*q_sq); //TODO: does the 2 need to be a 6?
                    for(int mu=0; mu<3; mu++)
                    {
                        noise_incr(i,j,k)(mu) = prefactor*get_rnd_gauss_fourier(i,j,k);
                    }
                    is_filled(i,j,k)=1.0;
                }
            }
        }
    }

    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                for(int mu=0; mu<3; mu++)
                {
                    xi_q(i,j,k)(mu) += (-dt)/tau*xi_q(i,j,k)(mu) + noise_incr(i,j,k)(mu);
                }
            }
        }
    }

}

void Generator::save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt)
{
    std::ofstream ofile;
    ofile.open(out_dir + "/noise_" + std::to_string(int(t)) + ".txt" );
    ofile << "nx: " << nx << std::endl;
    ofile << "ny: " << ny << std::endl;
    ofile << "nz: " << nz << std::endl;
    ofile << "dx: " << dx << std::endl;
    ofile << "time: " << (t*dt) << std::endl;
    ofile << "dt: " << dt << std::endl;
    ofile << "lambda: " << lambda << std::endl;
    ofile << "tau: " << tau << std::endl;
    ofile << "D: " << D << std::endl;
    ofile << "x y z xi_x xi_y xi_z" << std::endl;
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                ofile << i*dx << " " << j*dx << " " << k*dx << " " << theField(i,j,k)(0).real() << " " << theField(i,j,k)(1).real() << " " << theField(i,j,k)(2).real() << std::endl;
            }
        }
    }
}
#include "Generator.hpp"

Generator::Generator(ParamDict &theParams, gsl_rng *&the_rg)
{
    if(theParams.is_key("tau")) tau = std::stod(theParams.get_value("tau"));
    if(theParams.is_key("lambda")) lambda = std::stod(theParams.get_value("lambda"));
    if(theParams.is_key("D")) D = std::stod(theParams.get_value("D"));

    if(theParams.is_key("nx")) nx = std::stoi(theParams.get_value("nx"));
    if(theParams.is_key("ny")) ny = std::stoi(theParams.get_value("ny"));
    if(theParams.is_key("nz")) nz = std::stoi(theParams.get_value("nz"));

    if(theParams.is_key("dim")) dim = std::stoi(theParams.get_value("dim"));

    if(theParams.is_key("is_incompressible")) is_incompressible = std::stoi(theParams.get_value("is_incompressible"));

    for(int i=0; i<3; i++) spacing.push_back(1.0);
    if(theParams.is_key("dx")) spacing[0] = std::stod(theParams.get_value("dx"));
    if(theParams.is_key("dy")){
        spacing[1] = std::stod(theParams.get_value("dy"));
    }else if(theParams.is_key("dx")){
        spacing[1] = spacing[0];
    }
    if(theParams.is_key("dz")){
        spacing[2] = std::stod(theParams.get_value("dz"));
    } 
    else if(theParams.is_key("dx")){
        spacing[2] = spacing[0];
    }

    Lx = nx*spacing[0];
    Ly = ny*spacing[1];
    Lz = nz*spacing[2];

    //set RNG
    rg = the_rg;

    //Check compressibility is consistent with dim
    if (dim==1 && is_incompressible==1){
        std::cout << "WARNING: incompressibility cannot be enforced in 1d." << std::endl;
    }

    //Initialize noise differently for different no. of dimensions

    //3D
    if (dim==3) {
        xi_q = arma::field<arma::cx_vec>(nx, ny, nz);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    xi_q(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                }
            }
        }

        //White noise array in real space
        arma::cube real_white_noise_x = get_rnd_gauss_arr_3d(1.0);
        arma::cube real_white_noise_y = get_rnd_gauss_arr_3d(1.0);
        arma::cube real_white_noise_z = get_rnd_gauss_arr_3d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_cube fourier_white_noise_x = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_x);
        arma::cx_cube fourier_white_noise_y = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_y);
        arma::cx_cube fourier_white_noise_z = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_z);

        //Get spatial correlation function
        arma::cube cq = get_cq_3d();

        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    double prefactor = sqrt(D*Lx*Ly*Lz*cq(i,j,k));
                    xi_q(i,j,k)(0) = prefactor*fourier_white_noise_x(i,j,k);
                    xi_q(i,j,k)(1) = prefactor*fourier_white_noise_y(i,j,k);
                    xi_q(i,j,k)(2) = prefactor*fourier_white_noise_z(i,j,k);
                    if (is_incompressible==1){
                        //Get k vectors
                        arma::vec kVec(3, arma::fill::zeros);
                        kVec(0) = (2*M_PI/Lx)*(i-nx);
                        kVec(1) = (2*M_PI/Ly)*(j-ny);
                        kVec(2) = (2*M_PI/Lz)*(k-nz);
                        //Get projection operator
                        arma::mat projector = arma::mat(3,3,arma::fill::eye);
                        if (!(kVec(0)==0.0 && kVec(1)==0.0 && kVec(2)==0.0)){
                            arma::mat term2(3,3,arma::fill::zeros);
                            double denom = arma::dot(kVec, kVec);
                            for (int mu=0; mu<3; mu++){
                                for (int nu=0; nu<3; nu++){
                                    term2(mu,nu) = kVec(mu)*kVec(nu)/denom;
                                }
                            }
                            projector = projector - term2;
                        }
                        //Apply projection operator
                        xi_q(i,j,k) = projector*xi_q(i,j,k);
                    }
                }
            }
        }
    }

    //2D
    else if (dim==2) {
        xi_q = arma::field<arma::cx_vec>(nx, ny);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                xi_q(i,j) = arma::cx_vec(2, arma::fill::zeros);
            }
        }

        //White noise array in real space
        arma::mat real_white_noise_x = get_rnd_gauss_arr_2d(1.0);
        arma::mat real_white_noise_y = get_rnd_gauss_arr_2d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_mat fourier_white_noise_x = sqrt(nx*ny)*do_fourier_2d(real_white_noise_x);
        arma::cx_mat fourier_white_noise_y = sqrt(nx*ny)*do_fourier_2d(real_white_noise_y);

        //Get spatial correlation function
        arma::mat cq = get_cq_2d();

        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                double prefactor = sqrt(D*Lx*Ly*cq(i,j));
                xi_q(i,j)(0) = prefactor*fourier_white_noise_x(i,j);
                xi_q(i,j)(1) = prefactor*fourier_white_noise_y(i,j);
                if (is_incompressible==1){
                    //Get k vectors
                    arma::vec kVec(2, arma::fill::zeros);
                    kVec(0) = (2*M_PI/Lx)*(i-nx);
                    kVec(1) = (2*M_PI/Ly)*(j-ny);
                    //Get projection operator
                    arma::mat projector = arma::mat(2,2,arma::fill::eye);
                    if (!(kVec(0)==0.0 && kVec(1)==0.0)){
                        arma::mat term2(2,2,arma::fill::zeros);
                        double denom = arma::dot(kVec, kVec);
                        for (int mu=0; mu<2; mu++){
                            for (int nu=0; nu<2; nu++){
                                term2(mu,nu) = kVec(mu)*kVec(nu)/denom;
                            }
                        }
                        projector = projector - term2;
                    }
                    //Apply projection operator
                    xi_q(i,j) = projector*xi_q(i,j);
                }
            }
        }
    }

    //1D
    else if (dim==1) {
        xi_q = arma::field<arma::cx_vec>(nx);
        for(int i=0; i<nx; i++) {
            xi_q(i) = arma::cx_vec(1, arma::fill::zeros);
        }

        //White noise array in real space
        arma::vec real_white_noise_x = get_rnd_gauss_arr_1d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_vec fourier_white_noise_x = sqrt(nx)*do_fourier_1d(real_white_noise_x);

        //Get spatial correlation function
        arma::vec cq = get_cq_1d();

        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            double prefactor = sqrt(D*Lx*cq(i));
            xi_q(i)(0) = prefactor*fourier_white_noise_x(i);
        }
    }

    //Otherwise throw error
    else {
        std::cout << "Error: " << dim << "-dimensional noise field not supported." << std::endl;
        exit(-1);
    }
    
    std::cout << "Initialized " << dim << "-dimensional generator." << std::endl;
}

Generator::~Generator() {}

arma::field<arma::cx_vec> Generator::get_xi_r(int do_fft) {

    if(do_fft==1) return get_xi_r_fast();

    using namespace std::complex_literals;

    //3D
    if (dim==3) {
        arma::field<arma::cx_vec> xi_r(nx, ny, nz);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    xi_r(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                }
            }
        }

        //Not-so-fast Fourier transform
        for(int mu=0; mu<3; mu++) {
            for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                    for(int k=0; k<nz; k++) {
                        for(int q1=0; q1<nx; q1++) {
                            for(int q2=0; q2<ny; q2++) {
                                for(int q3=0; q3<nz; q3++) {
                                    std::complex<double> update = xi_q(q1,q2,q3)(mu)*std::exp(2*M_PI*1i*(
                                        (1.0*i*q1)/(1.0*nx) +
                                        (1.0*j*q2)/(1.0*ny) +
                                        (1.0*k*q3)/(1.0*nz)));
                                    //Check that imaginary part is zero
                                    if(update.imag()>1e-3){
                                        std::cout << "Warning: imaginary part of transform is not zero!" << std::endl;
                                        std::cout << "imaginary part:" << update.imag() << std::endl;
                                    }
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

    //2D
    else if (dim==2) {
        arma::field<arma::cx_vec> xi_r(nx, ny);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                xi_r(i,j) = arma::cx_vec(2, arma::fill::zeros);
            }
        }

        //Not-so-fast Fourier transform
        for(int mu=0; mu<2; mu++) {
            for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                    for(int q1=0; q1<nx; q1++) {
                        for(int q2=0; q2<ny; q2++) {
                            std::complex<double> update = xi_q(q1,q2)(mu)*std::exp(2*M_PI*1i*(
                                (1.0*i*q1)/(1.0*nx) +
                                (1.0*j*q2)/(1.0*ny)));
                            //Check that imaginary part is zero
                            /*
                            if(update.imag()>1e-3){
                                std::cout << "Warning: imaginary part of transform is not zero!" << std::endl;
                                std::cout << "imaginary part:" << update.imag() << std::endl;
                            }
                            */
                            xi_r(i,j)(mu) += update.real();
                        }
                    }
                    xi_r(i,j)(mu) *= 1.0/(Lx*Ly);
                }
            }
        }

        return xi_r;
    }

    //1D
    else if (dim==1) {
        arma::field<arma::cx_vec> xi_r(nx);
        for(int i=0; i<nx; i++) {
            xi_r(i) = arma::cx_vec(1, arma::fill::zeros);
        }

        //Not-so-fast Fourier transform
        for(int i=0; i<nx; i++) {
            for(int q1=0; q1<nx; q1++) {
                std::complex<double> update = xi_q(q1)(0)*std::exp(2*M_PI*1i*(1.0*i*q1)/(1.0*nx));
                //Check that imaginary part is zero
                if(update.imag()>1e-3){
                    std::cout << "Warning: imaginary part of transform is not zero!" << std::endl;
                    std::cout << "imaginary part:" << update.imag() << std::endl;
                }
                xi_r(i)(0) += update.real();
            }
            xi_r(i)(0) *= 1.0/(Lx);
        }

        return xi_r;
    }

    else {
        std::cout << "Error: " << dim << "-dimensional noise field not supported." << std::endl;
        exit(-1);
    }
}

arma::field<arma::cx_vec> Generator::get_xi_r_fast()
{
    using namespace std::complex_literals;

    //3D
    if (dim==3) {

        //Initialize noise field to zero
        arma::field<arma::cx_vec> xi_r(nx, ny, nz);
        arma::cx_cube xi_rx(nx, ny, nz, arma::fill::zeros);
        arma::cx_cube xi_ry(nx, ny, nz, arma::fill::zeros);
        arma::cx_cube xi_rz(nx, ny, nz, arma::fill::zeros);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    xi_rx(i,j,k) = xi_q(i,j,k)(0);
                    xi_ry(i,j,k) = xi_q(i,j,k)(1);
                    xi_rz(i,j,k) = xi_q(i,j,k)(2);
                    xi_r(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                }
            }
        }

        //Fast Fourier transform w/ FFTW
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

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
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

    //2D
    else if (dim==2) {

        //Initialize noise field to zero
        arma::field<arma::cx_vec> xi_r(nx, ny);
        arma::cx_mat xi_rx(nx, ny, arma::fill::zeros);
        arma::cx_mat xi_ry(nx, ny, arma::fill::zeros);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                xi_rx(i,j) = xi_q(i,j)(0);
                xi_ry(i,j) = xi_q(i,j)(1);
                xi_r(i,j) = arma::cx_vec(2, arma::fill::zeros);
            }
        }

        //Fast Fourier transform w/ FFTW
        fftw_complex *in_x, *in_y;
        fftw_plan px, py;

        //Might be faster to initialize the plan in the constructor since all these transforms have the same size
        //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
        in_x = reinterpret_cast<fftw_complex*>(xi_rx.memptr());
        in_y = reinterpret_cast<fftw_complex*>(xi_ry.memptr());
        px = fftw_plan_dft_2d(nx, ny, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);
        py = fftw_plan_dft_2d(nx, ny, in_y, in_y, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(px); 
        fftw_execute(py);

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                xi_r(i,j)(0) = xi_rx(i,j)/(Lx*Ly);
                xi_r(i,j)(1) = xi_ry(i,j)/(Lx*Ly);
            }
        }
        fftw_destroy_plan(px);
        fftw_destroy_plan(py);

        return xi_r;
    }

    //1D
    else if (dim==1) {

        //Initialize noise field to zero
        arma::field<arma::cx_vec> xi_r(nx);
        arma::cx_vec xi_rx(nx, arma::fill::zeros);
        for(int i=0; i<nx; i++) {
            xi_rx(i) = xi_q(i)(0);
            xi_r(i) = arma::cx_vec(1, arma::fill::zeros);
        }

        //Fast Fourier transform w/ FFTW
        fftw_complex *in_x;
        fftw_plan px;

        //Might be faster to initialize the plan in the constructor since all these transforms have the same size
        //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
        in_x = reinterpret_cast<fftw_complex*>(xi_rx.memptr());
        px = fftw_plan_dft_1d(nx, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(px); 

        for(int i=0; i<nx; i++) {
            xi_r(i)(0) = xi_rx(i)/(Lx);
        }
        fftw_destroy_plan(px);

        return xi_r;
    }

    else {
        std::cout << "Error: " << dim << "-dimensional noise field not supported." << std::endl;
        exit(-1);
    }
}

arma::cx_vec Generator::do_fourier_1d(arma::vec real_arr){
    
    arma::cx_vec fourier_arr(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        fourier_arr(i) = real_arr(i);
    }

    //Fast Fourier transform w/ FFTW
    fftw_complex *in_x;
    fftw_plan px;

    in_x = reinterpret_cast<fftw_complex*>(fourier_arr.memptr());
    px = fftw_plan_dft_1d(nx, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(px); 

    fourier_arr = fourier_arr/(nx);

    fftw_destroy_plan(px);

    return fourier_arr;
}

arma::cx_mat Generator::do_fourier_2d(arma::mat real_arr){
    
    arma::cx_mat fourier_arr(nx, ny, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            fourier_arr(i,j) = real_arr(i,j);
        }
    }

    //Fast Fourier transform w/ FFTW
    fftw_complex *in_x;
    fftw_plan px;

    in_x = reinterpret_cast<fftw_complex*>(fourier_arr.memptr());
    px = fftw_plan_dft_2d(nx, ny, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(px); 

    fourier_arr = fourier_arr/(nx*ny);

    fftw_destroy_plan(px);

    return fourier_arr;
}

arma::cx_cube Generator::do_fourier_3d(arma::cube real_arr){
    
    arma::cx_cube fourier_arr(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                fourier_arr(i,j,k) = real_arr(i,j,k);
            }
        }
    }

    //Fast Fourier transform w/ FFTW
    fftw_complex *in_x;
    fftw_plan px;

    in_x = reinterpret_cast<fftw_complex*>(fourier_arr.memptr());
    px = fftw_plan_dft_3d(nx, ny, nz, in_x, in_x, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(px); 

    fourier_arr = fourier_arr/(nx*ny*nz);

    fftw_destroy_plan(px);

    return fourier_arr;
}

std::complex<double> Generator::get_rnd_gauss_fourier_3D(int i, int j, int k)
{
    using namespace std::complex_literals;

    if( (i==0 && j==0 && k==0) ||
        (i==nx/2 && j==0 && k==0) ||
        (i==0 && j==ny/2 && k==0) ||
        (i==0 && j==0 && k==nz/2) ||
        (i==nx/2 && j==ny/2 && k==0) ||
        (i==nx/2 && j==0 && k==nz/2) ||
        (i==0 && j==ny/2 && k==nz/2) ||
        (i==nx/2 && j==ny/2 && k==nz/2) ) {

        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

std::complex<double> Generator::get_rnd_gauss_fourier_2D(int i, int j)
{
    using namespace std::complex_literals;

    if( (i==0 && j==0) ||
        (i==nx/2 && j==0) ||
        (i==0 && j==ny/2) ||
        (i==nx/2 && j==ny/2) ) {
        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

std::complex<double> Generator::get_rnd_gauss_fourier_1D(int i)
{
    using namespace std::complex_literals;

    if ( i==0 || i==nx/2 ) {
        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

arma::vec Generator::get_cq_1d(){

    //Get k vectors
    arma::vec kVec_x(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        if(i<=nx/2) kVec_x(i) = (2*M_PI/Lx)*i;
        else kVec_x(i) = (2*M_PI/Lx)*(i-nx);
    }

    //Compute normalized q-space correlation function
    arma::vec cq(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        double q_sq = kVec_x[i]*kVec_x[i];
        cq(i) = 2*lambda/(1+lambda*lambda*q_sq);
    }
    cq = cq/(arma::accu(cq)/(Lx)); //normalize such that f(r=0)=1
    if (fabs(arma::accu(cq)/(Lx)-1)>1e-10){
        std::cout << "Warning: c(q) is not normalized: " << arma::accu(cq)/(Lx) << std::endl;
    }

    return cq;
}

arma::mat Generator::get_cq_2d(){

    //Get k vectors
    arma::vec kVec_x(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        if(i<=nx/2) kVec_x(i) = (2*M_PI/Lx)*i;
        else kVec_x(i) = (2*M_PI/Lx)*(i-nx);
    }
    arma::vec kVec_y(ny, arma::fill::zeros);
    for(int i=0; i<ny; i++){
        if(i<=ny/2) kVec_y(i) = (2*M_PI/Ly)*i;
        else kVec_y(i) = (2*M_PI/Ly)*(i-ny);
    }

    //Compute normalized q-space correlation function
    arma::mat cq(nx, ny, arma::fill::zeros);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            double q_sq = kVec_x[i]*kVec_x[i] + kVec_y[j]*kVec_y[j];
            cq(i,j) = 2*M_PI*lambda*lambda*Lx*Ly/pow(1+lambda*lambda*q_sq,3.0/2.0);
        }
    }
    cq = cq/(arma::accu(cq)/(Lx*Ly)); //normalize such that f(r=0)=1
    if (fabs(arma::accu(cq)/(Lx*Ly)-1)>1e-10){
        std::cout << "Warning: c(q) is not normalized: " << arma::accu(cq)/(Lx*Ly) << std::endl;
    }

    return cq;
}

arma::cube Generator::get_cq_3d(){

    //Get k vectors
    arma::vec kVec_x(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        if(i<=nx/2) kVec_x(i) = (2*M_PI/Lx)*i;
        else kVec_x(i) = (2*M_PI/Lx)*(i-nx);
    }
    arma::vec kVec_y(ny, arma::fill::zeros);
    for(int i=0; i<ny; i++){
        if(i<=ny/2) kVec_y(i) = (2*M_PI/Ly)*i;
        else kVec_y(i) = (2*M_PI/Ly)*(i-ny);
    }
    arma::vec kVec_z(nz, arma::fill::zeros);
    for(int i=0; i<nz; i++){
        if(i<=nz/2) kVec_z(i) = (2*M_PI/Lz)*i;
        else kVec_z(i) = (2*M_PI/Lz)*(i-nz);
    }

    //Compute normalized q-space correlation function
    arma::cube cq(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                double q_sq = kVec_x[i]*kVec_x[i] + kVec_y[j]*kVec_y[j] + kVec_z[k]*kVec_z[k];
                cq(i,j,k) = 8*M_PI*lambda*lambda*lambda*Lx*Ly*Lz/((1+lambda*lambda*q_sq)*(1+lambda*lambda*q_sq));
            }
        }
    }
    cq = cq/(arma::accu(cq)/(Lx*Ly*Lz)); //normalize such that f(r=0)=1
    if (fabs(arma::accu(cq)/(Lx*Ly*Lz)-1)>1e-10){
        std::cout << "Warning: c(q) is not normalized: " << arma::accu(cq)/(Lx*Ly*Lz) << std::endl;
    }
    
    return cq;
}

arma::vec Generator::get_rnd_gauss_arr_1d(double var){

    arma::vec mymat(nx, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        mymat(i) = gsl_ran_gaussian(rg, var);
    }

    return mymat;
} 

arma::mat Generator::get_rnd_gauss_arr_2d(double var){

    arma::mat mymat(nx, ny, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            mymat(i,j) = gsl_ran_gaussian(rg, var);
        }
    }

    return mymat;
}

arma::cube Generator::get_rnd_gauss_arr_3d(double var){

    arma::cube mymat(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                mymat(i,j,k) = gsl_ran_gaussian(rg, var);
            }
        }
    }

    return mymat;
} 

void Generator::step(double dt)
{
    //3D
    if (dim==3) {
        //Initialize noise increment to zero
        arma::field<arma::cx_vec> noise_incr(nx, ny, nz);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    noise_incr(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
                }
            }
        }

        //White noise array in real space
        arma::cube real_white_noise_x = get_rnd_gauss_arr_3d(1.0);
        arma::cube real_white_noise_y = get_rnd_gauss_arr_3d(1.0);
        arma::cube real_white_noise_z = get_rnd_gauss_arr_3d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_cube fourier_white_noise_x = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_x);
        arma::cx_cube fourier_white_noise_y = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_y);
        arma::cx_cube fourier_white_noise_z = sqrt(nx*ny*nz)*do_fourier_3d(real_white_noise_z);

        //Get spatial correlation function
        arma::cube cq = get_cq_3d();

        //Compute noise increment
        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    double prefactor = sqrt(2*D*dt*Lx*Ly*Lz*cq(i,j,k)/tau);
                    noise_incr(i,j,k)(0) = prefactor*fourier_white_noise_x(i,j,k);
                    noise_incr(i,j,k)(1) = prefactor*fourier_white_noise_y(i,j,k);
                    noise_incr(i,j,k)(2) = prefactor*fourier_white_noise_z(i,j,k);
                    if (is_incompressible==1){
                        //Get k vectors
                        arma::vec kVec(3, arma::fill::zeros);
                        kVec(0) = (2*M_PI/Lx)*(i-nx);
                        kVec(1) = (2*M_PI/Ly)*(j-ny);
                        kVec(2) = (2*M_PI/Lz)*(k-nz);
                        //Get projection operator
                        arma::mat projector = arma::mat(3,3,arma::fill::eye);
                        if (!(kVec(0)==0.0 && kVec(1)==0.0 && kVec(2)==0.0)){
                            arma::mat term2(3,3,arma::fill::zeros);
                            double denom = arma::dot(kVec, kVec);
                            for (int mu=0; mu<3; mu++){
                                for (int nu=0; nu<3; nu++){
                                    term2(mu,nu) = kVec(mu)*kVec(nu)/denom;
                                }
                            }
                            projector = projector - term2;
                        }
                        //Apply projection operator
                        noise_incr(i,j,k) = projector*noise_incr(i,j,k);
                    }
                }
            }
        }

        //Update noise according to its OU dynamics
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    for(int mu=0; mu<3; mu++) {
                        xi_q(i,j,k)(mu) += (-dt)/tau*xi_q(i,j,k)(mu) + noise_incr(i,j,k)(mu);
                    }
                }
            }
        }
    }

    //2D
    else if (dim==2){
        //Initialize noise increment to zero
        arma::field<arma::cx_vec> noise_incr(nx, ny);
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                noise_incr(i,j) = arma::cx_vec(2, arma::fill::zeros);
            }
        }

        //White noise array in real space
        arma::mat real_white_noise_x = get_rnd_gauss_arr_2d(1.0);
        arma::mat real_white_noise_y = get_rnd_gauss_arr_2d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_mat fourier_white_noise_x = sqrt(nx*ny)*do_fourier_2d(real_white_noise_x);
        arma::cx_mat fourier_white_noise_y = sqrt(nx*ny)*do_fourier_2d(real_white_noise_y);

        //Get spatial correlation function
        arma::mat cq = get_cq_2d();

        //Compute noise increment
        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                double prefactor = sqrt(2*D*dt*Lx*Ly*cq(i,j)/tau);
                noise_incr(i,j)(0) = prefactor*fourier_white_noise_x(i,j);
                noise_incr(i,j)(1) = prefactor*fourier_white_noise_y(i,j);
                if (is_incompressible==1){
                    //Get k vectors
                    arma::vec kVec(2, arma::fill::zeros);
                    kVec(0) = (2*M_PI/Lx)*(i-nx);
                    kVec(1) = (2*M_PI/Ly)*(j-ny);
                    //Get projection operator
                    arma::mat projector = arma::mat(2,2,arma::fill::eye);
                    if (!(kVec(0)==0.0 && kVec(1)==0.0)){
                        arma::mat term2(2,2,arma::fill::zeros);
                        double denom = arma::dot(kVec, kVec);
                        for (int mu=0; mu<2; mu++){
                            for (int nu=0; nu<2; nu++){
                                term2(mu,nu) = kVec(mu)*kVec(nu)/denom;
                            }
                        }
                        projector = projector - term2;
                    }
                    //Apply projection operator
                    noise_incr(i,j) = projector*noise_incr(i,j);
                }
            }
        }

        //Update noise according to its OU dynamics
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int mu=0; mu<2; mu++) {
                    xi_q(i,j)(mu) += (-dt)/tau*xi_q(i,j)(mu) + noise_incr(i,j)(mu);
                    //xi_q(i,j)(mu) = exp(-dt/tau)*xi_q(i,j)(mu) + noise_incr(i,j)(mu); //''exact'' version
                }
            }
        }
    }

    //1D
    else if (dim==1){
        arma::cx_vec noise_incr(nx,arma::fill::zeros);

        //White noise array in real space
        arma::vec real_white_noise_x = get_rnd_gauss_arr_1d(1.0);

        //Take Fourier transform
        //Multiply by sqrt(N/2) to get unit variance in reciprocal space
        arma::cx_vec fourier_white_noise_x = sqrt(nx)*do_fourier_1d(real_white_noise_x);

        //Get spatial correlation function
        arma::vec cq = get_cq_1d();

        //Compute noise increment
        //Convolve white noise with correlation matrix
        for(int i=0; i<nx; i++) {
            double prefactor = sqrt(2*D*dt*Lx*cq(i)/tau);
            noise_incr(i) = prefactor*fourier_white_noise_x(i);
        }

        //Update noise according to its OU dynamics
        for(int i=0; i<nx; i++) {
            xi_q(i)(0) += (-dt)/tau*xi_q(i)(0) + noise_incr(i);
        }
    }

    else {
        std::cout << "Error: " << dim << "-dimensional noise field not supported." << std::endl;
        exit(-1);
    }
}

void Generator::open_h5(std::string out_dir) {
    //Create an empty h5 file for storing noise trajectory
    using namespace HighFive;
    fs::create_directories(out_dir);
    std::string name = out_dir + "/noise_traj.h5";
    if(fs::exists(name)) {
        std::cout << "Warning: file already exists. Overwriting..." << std::endl; 
        fs::remove(name);
    }
    File file(name, File::ReadWrite | File::Create | File::Truncate);

    //Data groups
    Group parameters = file.createGroup("/parameters");
    Group grid = file.createGroup("/grid");
    Group noise = file.createGroup("/noise");

    //Write parameters (which don't change over trajectory)
    if (dim==3) {
        std::vector<int> dims{nx, ny, nz};
        DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
        dimensions.write(dims);
        DataSet spacing_data = file.createDataSet<double>("/grid/spacing", DataSpace::From(spacing));
        spacing_data.write(spacing);
    }
    else if (dim==2) {
        std::vector<int> dims{nx, ny};
        DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
        dimensions.write(dims);
        DataSet spacing_data = file.createDataSet<double>("/grid/spacing", DataSpace::From(spacing));
        spacing_data.write(spacing);
    }
    else if (dim==1) {
        std::vector<int> dims{nx};
        DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
        dimensions.write(dims);
        DataSet spacing_data = file.createDataSet<double>("/grid/spacing", DataSpace::From(spacing));
        spacing_data.write(spacing);    
    }


    DataSet lambda_value = file.createDataSet<double>("/parameters/lambda", DataSpace::From(lambda));
    lambda_value.write(lambda);

    DataSet tau_value = file.createDataSet<double>("/parameters/tau", DataSpace::From(tau));
    tau_value.write(tau);

    DataSet D_value = file.createDataSet<double>("/parameters/D", DataSpace::From(D));
    D_value.write(D);
}

void Generator::save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt)
{
    double t_curr = t*dt;
    using namespace HighFive;

    if (dim==3) {
        try {
            std::string name = out_dir + "/noise_traj.h5";
            if(!fs::exists(name)) {
                std::cout << "Error: file does not exist!" << std::endl; 
                exit(0);
            }
            File file(name, File::ReadWrite);

            //Create necessary data structures for storage
            //( HighFive may not understand armadillo)
            std::vector<std::vector<std::vector<std::vector<double>>>> noise_val_x(1, std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny,std::vector<double>(nz,0.0))));
            std::vector<std::vector<std::vector<std::vector<double>>>> noise_val_y(1, std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny,std::vector<double>(nz,0.0))));
            std::vector<std::vector<std::vector<std::vector<double>>>> noise_val_z(1, std::vector<std::vector<std::vector<double>>>(nx, std::vector<std::vector<double>>(ny,std::vector<double>(nz,0.0))));

            //Fill in noise values
            for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                    for(int k=0; k<nz; k++) {
                        noise_val_x[0][i][j][k] = theField(i,j,k)(0).real();
                        noise_val_y[0][i][j][k] = theField(i,j,k)(1).real();
                        noise_val_z[0][i][j][k] = theField(i,j,k)(2).real();
                    }
                }
            }

            //Allocate memory for data
            DataSpace part_t_space = DataSpace({1},{DataSpace::UNLIMITED});
            DataSetCreateProps props_time;
            props_time.add(Chunking(std::vector<hsize_t>{1}));

            DataSpace part_val_space = DataSpace({1, nx, ny, nz},{DataSpace::UNLIMITED, nx, ny, nz});
            DataSetCreateProps props_val;
            props_val.add(Chunking(std::vector<hsize_t>{1,nx,ny,nz}));

            //Create datasets if it is the first time writing to this trajectory
            if(!file.exist("/noise/time")) {

                //Create time data
                DataSet time = file.createDataSet<double>("/noise/time", part_t_space, props_time);

                time.write(t_curr);

                //Create noise data
                DataSet value_x = file.createDataSet<double>("/noise/value/x", part_val_space, props_val);
                DataSet value_y = file.createDataSet<double>("/noise/value/y", part_val_space, props_val);
                DataSet value_z = file.createDataSet<double>("/noise/value/z", part_val_space, props_val);

                value_x.select({0,0,0,0},{1,nx,ny,nz}).write(noise_val_x);
                value_y.select({0,0,0,0},{1,nx,ny,nz}).write(noise_val_y);
                value_z.select({0,0,0,0},{1,nx,ny,nz}).write(noise_val_z);
            }
            //Otherwise, append to the trajectory
            else {
                
                //Update time
                DataSet time = file.getDataSet("/noise/time");
                std::vector<long unsigned int> time_dim = time.getDimensions();
                std::vector<long unsigned int> time_dim_old = time_dim;
                time_dim[0] += 1;
                time.resize(time_dim);
                time.select(time_dim_old,{1}).write(t_curr);

                //Update noise
                DataSet value_x = file.getDataSet("/noise/value/x");
                std::vector<long unsigned int> value_x_dim = value_x.getDimensions();
                std::vector<long unsigned int> value_x_dim_old = value_x_dim;
                value_x_dim[0] += 1;
                value_x.resize(value_x_dim);
                value_x.select({value_x_dim_old[0],0,0,0},{1,nx,ny,nz}).write(noise_val_x);

                DataSet value_y = file.getDataSet("/noise/value/y");
                std::vector<long unsigned int> value_y_dim = value_y.getDimensions();
                std::vector<long unsigned int> value_y_dim_old = value_y_dim;
                value_y_dim[0] += 1;
                value_y.resize(value_y_dim);
                value_y.select({value_y_dim_old[0],0,0,0},{1,nx,ny,nz}).write(noise_val_y);

                DataSet value_z = file.getDataSet("/noise/value/z");
                std::vector<long unsigned int> value_z_dim = value_z.getDimensions();
                std::vector<long unsigned int> value_z_dim_old = value_z_dim;
                value_z_dim[0] += 1;
                value_z.resize(value_z_dim);
                value_z.select({value_z_dim_old[0],0,0,0},{1,nx,ny,nz}).write(noise_val_z);
            }
        }
        catch(Exception& err) {
            std::cerr << err.what() << std::endl;
        }
    }

    else if (dim==2) {
        try {
            std::string name = out_dir + "/noise_traj.h5";
            if(!fs::exists(name)) {
                std::cout << "Error: file does not exist!" << std::endl; 
                exit(0);
            }
            File file(name, File::ReadWrite);

            //Create necessary data structures for storage
            //( HighFive may not understand armadillo)
            std::vector<std::vector<std::vector<double>>> noise_val_x(1, std::vector<std::vector<double>>(nx, std::vector<double>(ny,0.0)));
            std::vector<std::vector<std::vector<double>>> noise_val_y(1, std::vector<std::vector<double>>(nx, std::vector<double>(ny,0.0)));

            //Fill in noise values
            for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                    noise_val_x[0][i][j] = theField(i,j)(0).real();
                    noise_val_y[0][i][j] = theField(i,j)(1).real();
                }
            }

            //Allocate memory for data
            DataSpace part_t_space = DataSpace({1},{DataSpace::UNLIMITED});
            DataSetCreateProps props_time;
            props_time.add(Chunking(std::vector<hsize_t>{1}));

            DataSpace part_val_space = DataSpace({1, nx, ny},{DataSpace::UNLIMITED, nx, ny});
            DataSetCreateProps props_val;
            props_val.add(Chunking(std::vector<hsize_t>{1,nx,ny}));

            //Create datasets if it is the first time writing to this trajectory
            if(!file.exist("/noise/time")) {

                //Create time data
                DataSet time = file.createDataSet<double>("/noise/time", part_t_space, props_time);

                time.write(t_curr);

                //Create noise data
                DataSet value_x = file.createDataSet<double>("/noise/value/x", part_val_space, props_val);
                DataSet value_y = file.createDataSet<double>("/noise/value/y", part_val_space, props_val);

                value_x.select({0,0,0},{1,nx,ny}).write(noise_val_x);
                value_y.select({0,0,0},{1,nx,ny}).write(noise_val_y);
            }
            //Otherwise, append to the trajectory
            else {
                
                //Update time
                DataSet time = file.getDataSet("/noise/time");
                std::vector<long unsigned int> time_dim = time.getDimensions();
                std::vector<long unsigned int> time_dim_old = time_dim;
                time_dim[0] += 1;
                time.resize(time_dim);
                time.select(time_dim_old,{1}).write(t_curr);

                //Update noise
                DataSet value_x = file.getDataSet("/noise/value/x");
                std::vector<long unsigned int> value_x_dim = value_x.getDimensions();
                std::vector<long unsigned int> value_x_dim_old = value_x_dim;
                value_x_dim[0] += 1;
                value_x.resize(value_x_dim);
                value_x.select({value_x_dim_old[0],0,0},{1,nx,ny}).write(noise_val_x);

                DataSet value_y = file.getDataSet("/noise/value/y");
                std::vector<long unsigned int> value_y_dim = value_y.getDimensions();
                std::vector<long unsigned int> value_y_dim_old = value_y_dim;
                value_y_dim[0] += 1;
                value_y.resize(value_y_dim);
                value_y.select({value_y_dim_old[0],0,0},{1,nx,ny}).write(noise_val_y);
            }
        }
        catch(Exception& err) {
            std::cerr << err.what() << std::endl;
        }
    }

    else if (dim==1) {
        try {
            std::string name = out_dir + "/noise_traj.h5";
            if(!fs::exists(name)) {
                std::cout << "Error: file does not exist!" << std::endl; 
                exit(0);
            }
            File file(name, File::ReadWrite);

            //Create necessary data structures for storage
            //( HighFive may not understand armadillo)
            std::vector<std::vector<double>> noise_val_x(1, std::vector<double>(nx, 0.0));

            //Fill in noise values
            for(int i=0; i<nx; i++) {
                noise_val_x[0][i] = theField(i)(0).real();
            }

            //Allocate memory for data
            DataSpace part_t_space = DataSpace({1},{DataSpace::UNLIMITED});
            DataSetCreateProps props_time;
            props_time.add(Chunking(std::vector<hsize_t>{1}));

            DataSpace part_val_space = DataSpace({1, nx},{DataSpace::UNLIMITED, nx});
            DataSetCreateProps props_val;
            props_val.add(Chunking(std::vector<hsize_t>{1,nx}));

            //Create datasets if it is the first time writing to this trajectory
            if(!file.exist("/noise/time")) {

                //Create time data
                DataSet time = file.createDataSet<double>("/noise/time", part_t_space, props_time);

                time.write(t_curr);

                //Create noise data
                DataSet value_x = file.createDataSet<double>("/noise/value/x", part_val_space, props_val);

                value_x.select({0,0},{1,nx}).write(noise_val_x);
            }
            //Otherwise, append to the trajectory
            else {
                
                //Update time
                DataSet time = file.getDataSet("/noise/time");
                std::vector<long unsigned int> time_dim = time.getDimensions();
                std::vector<long unsigned int> time_dim_old = time_dim;
                time_dim[0] += 1;
                time.resize(time_dim);
                time.select(time_dim_old,{1}).write(t_curr);

                //Update noise
                DataSet value_x = file.getDataSet("/noise/value/x");
                std::vector<long unsigned int> value_x_dim = value_x.getDimensions();
                std::vector<long unsigned int> value_x_dim_old = value_x_dim;
                value_x_dim[0] += 1;
                value_x.resize(value_x_dim);
                value_x.select({value_x_dim_old[0],0},{1,nx}).write(noise_val_x);
            }
        }
        catch(Exception& err) {
            std::cerr << err.what() << std::endl;
        }
    }
}

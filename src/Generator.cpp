#include "Generator.hpp"

/*****************************/
/* 3D Generator */
/*****************************/

Generator3D::Generator3D(ParamDict &theParams, gsl_rng *&the_rg)
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

    //set RNG
    rg = the_rg;

    xi_q = arma::field<arma::cx_vec>(nx, ny, nz);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                xi_q(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
            }
        }
    }

    //Initialize q field to ensure stationarity
    //create a cube to track which field values have been filled
    arma::cube is_filled(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                if(is_filled((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)==1) {                    
                    for(int mu=0; mu<3; mu++) {
                        xi_q(i,j,k)(mu) = std::conj(xi_q((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)(mu));
                    }
                    is_filled(i,j,k)=1.0;
                }
                else {
                    double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
                    double cq_sqrt = 1/(1+lambda*lambda*q_sq);
                    double prefactor = sqrt(D/tau)*sqrt(Lx*Ly*Lz)*cq_sqrt;
                    for(int mu=0; mu<3; mu++){
                        xi_q(i,j,k)(mu) = prefactor*get_rnd_gauss_fourier(i,j,k);
                    }
                    is_filled(i,j,k)=1.0;
                }
            }
        }
    }
    std::cout << "Initialized generator." << std::endl;
}

Generator3D::~Generator3D() {}

arma::field<arma::cx_vec> Generator3D::get_xi_r(int do_fft)
{

    if(do_fft==1) return get_xi_r_fast();

    using namespace std::complex_literals;

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

arma::field<arma::cx_vec> Generator3D::get_xi_r_fast()
{
    using namespace std::complex_literals;

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

std::complex<double> Generator3D::get_rnd_gauss_fourier(int i, int j, int k)
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

void Generator3D::step(double dt)
{
    //Initialize noise increment to zero
    arma::field<arma::cx_vec> noise_incr(nx, ny, nz);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                noise_incr(i,j,k) = arma::cx_vec(3, arma::fill::zeros);
            }
        }
    }

    //Compute noise increment
    arma::cube is_filled(nx, ny, nz, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                //Check for pairs of values that are complex conjugates by symmetry
                if(is_filled((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)==1) {
                    for(int mu=0; mu<3; mu++) {
                        noise_incr(i,j,k)(mu) = std::conj(noise_incr((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)(mu));
                    }
                    is_filled(i,j,k)=1.0;
                }
                else
                {
                    double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
                    //Assume that exponential is normalized, so prefactor is (8*M_PI*lambda*lambda*lambda)
                    double cq_sqrt = 1/(1+lambda*lambda*q_sq);
                    double prefactor = sqrt(2*D*dt*Lx*Ly*Lz)*cq_sqrt/tau;
                    for(int mu=0; mu<3; mu++) {
                        noise_incr(i,j,k)(mu) = prefactor*get_rnd_gauss_fourier(i,j,k);
                    }
                    is_filled(i,j,k)=1.0;
                }
            }
        }
    }

    //Update noise according to its OU dynamics
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

void Generator3D::open_h5(std::string out_dir) {
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
    std::vector<int> dims{nx, ny, nz};
    DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
    dimensions.write(dims);

    DataSet spacing = file.createDataSet<double>("/grid/spacing", DataSpace::From(dx));
    spacing.write(dx);

    DataSet lambda_value = file.createDataSet<double>("/parameters/lambda", DataSpace::From(lambda));
    lambda_value.write(lambda);

    DataSet tau_value = file.createDataSet<double>("/parameters/tau", DataSpace::From(tau));
    tau_value.write(tau);

    DataSet D_value = file.createDataSet<double>("/parameters/D", DataSpace::From(D));
    D_value.write(D);
}

void Generator3D::save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt)
{
    double t_curr = t*dt;
    using namespace HighFive;
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

/*****************************/
/* 2D Generator */
/*****************************/

Generator2D::Generator2D(ParamDict &theParams, gsl_rng *&the_rg)
{
    if(theParams.is_key("dx")) dx = std::stod(theParams.get_value("dx"));
    if(theParams.is_key("tau")) tau = std::stod(theParams.get_value("tau"));
    if(theParams.is_key("lambda")) lambda = std::stod(theParams.get_value("lambda"));
    if(theParams.is_key("D")) D = std::stod(theParams.get_value("D"));

    if(theParams.is_key("nx")) nx = std::stoi(theParams.get_value("nx"));
    if(theParams.is_key("ny")) ny = std::stoi(theParams.get_value("ny"));

    Lx = nx*dx;
    Ly = ny*dx;

    //set RNG
    rg = the_rg;

    xi_q = arma::field<arma::cx_vec>(nx, ny);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            xi_q(i,j) = arma::cx_vec(2, arma::fill::zeros);
        }
    }

    //Initialize q field to ensure stationarity
    //create a cube to track which field values have been filled
    arma::mat is_filled(nx, ny, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            if(is_filled((nx-i)%nx,(ny-j)%ny)==1) {                    
                for(int mu=0; mu<2; mu++) {
                    xi_q(i,j)(mu) = std::conj(xi_q((nx-i)%nx,(ny-j)%ny)(mu));
                }
                is_filled(i,j)=1.0;
            }
            else {
                double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly));
                //C(q) = 1/(1+q^2+lambda^2)^(3/2)
                double cq_sqrt = sqrt(1/((1+lambda*lambda*q_sq)*sqrt((1+lambda*lambda*q_sq))));
                double prefactor = sqrt(D/tau)*sqrt(Lx*Ly)*cq_sqrt;
                for(int mu=0; mu<2; mu++){
                    xi_q(i,j)(mu) = prefactor*get_rnd_gauss_fourier(i,j);
                }
                is_filled(i,j)=1.0;
            }
        }
    }
    std::cout << "Initialized generator." << std::endl;
}

Generator2D::~Generator2D() {}

arma::field<arma::cx_vec> Generator2D::get_xi_r(int do_fft)
{

    if(do_fft==1) return get_xi_r_fast();

    using namespace std::complex_literals;

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
                        if(update.imag()>1e-3){
                            std::cout << "Warning: imaginary part of transform is not zero!" << std::endl;
                            std::cout << "imaginary part:" << update.imag() << std::endl;
                        }
                        xi_r(i,j)(mu) += update.real();
                    }
                }
                xi_r(i,j)(mu) *= 1.0/(Lx*Ly);
            }
        }
    }
    return xi_r;
}

arma::field<arma::cx_vec> Generator2D::get_xi_r_fast()
{
    using namespace std::complex_literals;

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

std::complex<double> Generator2D::get_rnd_gauss_fourier(int i, int j)
{
    using namespace std::complex_literals;

    if( (i==0 && j==0) ||
        (i==nx/2 && j==0) ||
        (i==0 && j==ny/2) ||
        (i==nx/2 && j==ny/2) )
    {
        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else
    {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

void Generator2D::step(double dt)
{
    //Initialize noise increment to zero
    arma::field<arma::cx_vec> noise_incr(nx, ny);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            noise_incr(i,j) = arma::cx_vec(2, arma::fill::zeros);
        }
    }

    //Compute noise increment
    arma::mat is_filled(nx, ny, arma::fill::zeros);
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            //Check for pairs of values that are complex conjugates by symmetry
            if(is_filled((nx-i)%nx,(ny-j)%ny)==1) {
                for(int mu=0; mu<2; mu++) {
                    noise_incr(i,j)(mu) = std::conj(noise_incr((nx-i)%nx,(ny-j)%ny)(mu));
                }
                is_filled(i,j)=1.0;
            }
            else {
                double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly));
                double cq_sqrt = sqrt(1/((1+lambda*lambda*q_sq)*sqrt((1+lambda*lambda*q_sq))));
                double prefactor = sqrt(2*D*dt*Lx*Ly)*cq_sqrt/tau;
                for(int mu=0; mu<2; mu++) {
                    noise_incr(i,j)(mu) = prefactor*get_rnd_gauss_fourier(i,j);
                }
                is_filled(i,j)=1.0;
            }
        }
    }

    //Update noise according to its OU dynamics
    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int mu=0; mu<2; mu++) {
                xi_q(i,j)(mu) += (-dt)/tau*xi_q(i,j)(mu) + noise_incr(i,j)(mu);
            }
        }
    }

}

void Generator2D::open_h5(std::string out_dir) {
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
    std::vector<int> dims{nx, ny};
    DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
    dimensions.write(dims);

    DataSet spacing = file.createDataSet<double>("/grid/spacing", DataSpace::From(dx));
    spacing.write(dx);

    DataSet lambda_value = file.createDataSet<double>("/parameters/lambda", DataSpace::From(lambda));
    lambda_value.write(lambda);

    DataSet tau_value = file.createDataSet<double>("/parameters/tau", DataSpace::From(tau));
    tau_value.write(tau);

    DataSet D_value = file.createDataSet<double>("/parameters/D", DataSpace::From(D));
    D_value.write(D);
}

void Generator2D::save_field(arma::field<arma::cx_vec> &theField, std::string out_dir, double t, double dt)
{
    double t_curr = t*dt;
    using namespace HighFive;
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

/*****************************/
/* 1D Generator */
/*****************************/

Generator1D::Generator1D(ParamDict &theParams, gsl_rng *&the_rg)
{
    if(theParams.is_key("dx")) dx = std::stod(theParams.get_value("dx"));
    if(theParams.is_key("tau")) tau = std::stod(theParams.get_value("tau"));
    if(theParams.is_key("lambda")) lambda = std::stod(theParams.get_value("lambda"));
    if(theParams.is_key("D")) D = std::stod(theParams.get_value("D"));

    if(theParams.is_key("nx")) nx = std::stoi(theParams.get_value("nx"));

    Lx = nx*dx;

    xi_q = arma::cx_vec(nx,arma::fill::zeros);

    //set RNG
    rg = the_rg;

    //Initialize q field to ensure stationarity
    for(int i=0; i<nx/2+1; i++)
    {
        double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx));
        double prefactor = sqrt(D/tau)*sqrt(Lx)/sqrt(1+lambda*lambda*q_sq);

        xi_q(i) = prefactor*get_rnd_gauss_fourier(i);
        if(i>0) xi_q(nx-i) = std::conj(xi_q(i));
    }
}

Generator1D::~Generator1D() {}

arma::cx_vec Generator1D::get_xi_r(int do_fft)
{
    if(do_fft==1) return get_xi_r_fast();

    using namespace std::complex_literals;

    arma::cx_vec xi_r(nx,arma::fill::zeros);

    //Not-so-fast Fourier transform
    for(int i=0; i<nx; i++)
    {
        for(int q1=0; q1<nx; q1++)
        {
            std::complex<double> update = xi_q(q1)*std::exp(2*M_PI*1i*(1.0*i*q1/nx));
            xi_r(i) += update;
        }
        xi_r(i) *= (1.0/nx);
    }
    //std::cout << xi_r << std::endl;
    return xi_r;
}

arma::cx_vec Generator1D::get_xi_r_fast()
{
    using namespace std::complex_literals;

    arma::cx_vec xi_r(nx,arma::fill::zeros);
    for(int i=0; i<nx; i++) xi_r[i] = xi_q[i];

    //Fast Fourier transform w/ FFTW

    //std::cout << "doing fft" << std::endl;

    fftw_complex *in;
    fftw_plan p;

    //Might be faster to initialize the plan in the constructor since all these transforms have the same size
    //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    in = reinterpret_cast<fftw_complex*>(xi_r.memptr());
    p = fftw_plan_dft_1d(nx, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p); 

    for(int i=0; i<nx; i++)
    {
        xi_r(i) *= (1.0/nx);
    }
    //std::cout << xi_r << std::endl;
    fftw_destroy_plan(p);
    return xi_r;
}

std::complex<double> Generator1D::get_rnd_gauss_fourier(int i)
{
    using namespace std::complex_literals;

    if( i==0 || i==nx/2 )
    {
        return gsl_ran_gaussian(rg, 1.0) + 0i;
    }
    else
    {
        return gsl_ran_gaussian(rg, 1.0)/sqrt(2.0) + 1i*gsl_ran_gaussian(rg, 1.0)/sqrt(2.0);
    }
}

void Generator1D::step(double dt)
{
    arma::cx_vec noise_incr(nx);

    for(int i=0; i<nx/2+1; i++)
    {
        double q_sq = 4*M_PI*M_PI*(i*i/(Lx*Lx));
        double prefactor = sqrt(2*D*dt)/tau*sqrt(Lx)/sqrt(1+lambda*lambda*q_sq);
        noise_incr(i) = prefactor*get_rnd_gauss_fourier(i);
        if(i>0) noise_incr(nx-i) = std::conj(noise_incr(i));
    }

    for(int i=0; i<nx; i++)
    {
        xi_q(i) += (-dt)/tau*xi_q(i) + noise_incr(i);
    }

}

void Generator1D::open_h5(std::string out_dir) {
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
    std::vector<int> dims{nx};
    DataSet dimensions = file.createDataSet<int>("/grid/dimensions", DataSpace::From(dims));
    dimensions.write(dims);

    DataSet spacing = file.createDataSet<double>("/grid/spacing", DataSpace::From(dx));
    spacing.write(dx);

    DataSet lambda_value = file.createDataSet<double>("/parameters/lambda", DataSpace::From(lambda));
    lambda_value.write(lambda);

    DataSet tau_value = file.createDataSet<double>("/parameters/tau", DataSpace::From(tau));
    tau_value.write(tau);

    DataSet D_value = file.createDataSet<double>("/parameters/D", DataSpace::From(D));
    D_value.write(D);
}

void Generator1D::save_field(arma::cx_vec &theField, std::string out_dir, double t, double dt)
{
    double t_curr = t*dt;
    using namespace HighFive;
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
            noise_val_x[0][i] = theField(i).real();
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
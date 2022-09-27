/* Compute spatial and time correlation functions for active noise */

#include <cstdlib>
#include <iostream>
#include <filesystem>

#include "analysis.hpp"
#include "ParamDict.hpp"

int main(int argc, char *argv[])
{
    std::string conf_file = "analysis.conf";
    if (argc>1) conf_file = argv[1];

    ParamDict myParams;
    myParams.read_params(conf_file);

    std::string input_dir = "data";
    std::string output_dir = "example";
    double rmax=10.0;
    double delta_x = 1.0;
    double tmax=10.0;
    int nframes = 100;
    int nchunks = 1;
    int nframes_per_chunk;

    int nx;
    int ny;
    int nz;
    double dx;
    double tcurr;
    double dt;
    double lambda;
    double tau;
    double D;

    double Lx;
    double Ly;
    double Lz;

    if(myParams.is_key("input_dir")) input_dir = myParams.get_value("input_dir");
    if(myParams.is_key("output_dir")) output_dir = myParams.get_value("output_dir");
    if(myParams.is_key("rmax")) rmax = std::stod(myParams.get_value("rmax"));
    if(myParams.is_key("delta_x")) delta_x = std::stod(myParams.get_value("delta_x"));
    if(myParams.is_key("tmax")) tmax = std::stod(myParams.get_value("tmax"));
    if(myParams.is_key("nframes")) nframes = std::stoi(myParams.get_value("nframes"));
    if(myParams.is_key("nchunks")) nchunks = std::stoi(myParams.get_value("nchunks"));

    nframes_per_chunk = nframes/nchunks;
    std::cout << "nframes per chunk: " << nframes_per_chunk << std::endl;

    //Load parameters from configurations
    std::ifstream infile(input_dir + "/noise_0.txt");
    //TODO: yell if this file doesn't exist
    if (infile.is_open()) {
        std::string line;
        for(int i=0; i<9; i++) //TODO: change this so that it doesn't depend on how many params you have
        {
            std::getline(infile,line);
            int index = line.find(":");
            std::string key = line.substr(0,index);
            int index2 = line.find(" ");
            std::string value = line.substr(index2+1);
            if(key=="nx") nx = std::stoi(value);
            if(key=="ny") ny = std::stoi(value);
            if(key=="nz") nz = std::stoi(value);
            if(key=="dx") dx = std::stod(value);
            //if(key=="time") tmax = std::stod(value);
            if(key=="dt") dt = std::stod(value);
            if(key=="lambda") lambda = std::stod(value);
            if(key=="tau") tau = std::stod(value);
            if(key=="D") D = std::stod(value);
        }
        infile.close();
    } 

    Lx = nx*dx;
    Ly = ny*dx;
    Lz = nz*dx;
    arma::field<arma::vec> r(nx,ny,nz);
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                r(i,j,k) = arma::vec(3);
                r(i,j,k)(0) = i*dx;
                r(i,j,k)(1) = j*dx;
                r(i,j,k)(2) = k*dx;
            }
        }
    }

    std::cout << "tau: " << tau << std::endl;
    std::cout << "lambda: " << lambda << std::endl;
    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    std::cout << "nz: " << nz << std::endl;
    std::cout << "tmax: " << tmax << std::endl;
    std::cout << "dt: " << dt << std::endl;

    //Divide trajectory into chunks
    //(memory issue with loading entire large trajectories)
    std::cout << "Loading noise trajectory..." << std::endl;
    for(int c=0; c<nchunks; c++)
    {
        std::cout << "Analyzing chunk " << c << std::endl;
        std::vector<arma::field<arma::vec>> xi(nframes_per_chunk);
        for(int t=0; t<nframes_per_chunk; t++)
        {
            xi[t] = arma::field<arma::vec>(nx,ny,nz);
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                        xi[t](i,j,k) = arma::vec(3,arma::fill::zeros);
                    }
                }
            }
        }

        for(int t=0; t<nframes_per_chunk; t++)
        {
            if (t%1000==0) std::cout << "frame " << t+c*nframes_per_chunk << std::endl;
            //std::cout << "noise " << t+c*nframes_per_chunk << std::endl;
            std::ifstream infile2(input_dir + "/noise_" + std::to_string(t+c*nframes_per_chunk) + ".txt");
            int start_reading = 0;
            int count = 0;
            if (infile2.is_open()) {
                std::string line;
                while (std::getline(infile2, line)) 
                {
                    if(start_reading==1)
                    {
                        istringstream my_stream(line);
                        double x;
                        double y;
                        double z;
                        my_stream >> x;
                        my_stream >> y;
                        my_stream >> z;
                        my_stream >> xi[t](int(x/dx),int(y/dx),int(z/dx))(0);
                        my_stream >> xi[t](int(x/dx),int(y/dx),int(z/dx))(1);
                        my_stream >> xi[t](int(x/dx),int(y/dx),int(z/dx))(2);
                        count ++;                    
                    }

                    if(line.find(":") == std::string::npos && start_reading==0){
                        start_reading=1;
                        //std::cout << line << std::endl;
                    }
                }
                infile2.close();
            }
        }
        std::cout << xi[0].n_rows << " " << xi[0].n_cols << " " << xi[0].n_slices << std::endl;

        //Create vector of indices at which separation to compute
        //the spatial correlation function.
        arma::vec r1 = {0,0,0};
        std::vector<std::array<int, 3>> spat_indices = std::vector<std::array<int, 3>>();
        for(int i=0; i<nx; i++)
        {
            for(int j=0; j<ny; j++)
            {
                for(int k=0; k<nz; k++)
                {
                    arma::vec r2 = {dx*i, dx*j, dx*k};
                    double dist = get_min_image_dist(r1, r2, Lx, Ly, Lz);
                    if(dist<=rmax)
                    {
                        arma::vec disp = get_min_image_disp(r1, r2, Lx, Ly, Lz);
                        std::array<int, 3> arr;
                        for(int mu=0; mu<3; mu++) arr[mu] = (int)disp(mu)/dx;
                        spat_indices.push_back(arr);
                    }
                }
            }
        }
        //std::cout << spat_indices << std::endl;

        //Get correlation functions
        
        compute_spatial_corr(rmax, delta_x, Lx, Ly, Lz, xi, spat_indices, "_chunk=" + std::to_string(c), output_dir);
        compute_time_corr(tmax, dt, xi, "_chunk=" + std::to_string(c), output_dir);
    }
    
    return 0;
}

void compute_spatial_corr(double rmax, double delta_x, double Lx, double Ly, double Lz, std::vector<arma::field<arma::vec>> &xi, std::vector<std::array<int, 3>> &indices, std::string name_add, std::string out_dir)
{
    int nframes = int(xi.size());
    std::cout << "nframes: " << nframes << std::endl;
    int nx = xi[0].n_rows;
    int ny = xi[0].n_cols;
    int nz = xi[0].n_slices;
    double dx = Lx/nx;
    arma::vec r(int(rmax/delta_x)+1);
    for(int i=0; i<int(rmax/delta_x)+1; i++) r(i) = delta_x*i;
    arma::vec rhist(int(rmax/delta_x)+1,arma::fill::zeros);
    arma::vec c_r(int(rmax/delta_x)+1, arma::fill::zeros);
    int indices_size = int(indices.size());
    std::cout << "indices_size: " << indices_size << std::endl;
    arma::vec c_r_3d(indices_size, arma::fill::zeros);
    arma::vec rhist_3d(indices_size, arma::fill::zeros);

    std::cout << "Computing spatial correlation function..." << std::endl;
    for(int t=0; t<nframes; t++)
    {
        if (t%100==0) std::cout << "frame " << t << std::endl;
        for(int i1=0; i1<nx; i1++)
        {
            for(int j1=0; j1<ny; j1++)
            {
                for(int k1=0; k1<nz; k1++)
                {
                    //arma::vec r1 = {dx*i1, dx*j1, dx*k1};
                    for(int ind=0; ind<indices_size; ind++)
                    {
                        int i2 = i1+indices[ind][0];  
                        int j2 = j1+indices[ind][1];
                        int k2 = k1+indices[ind][2];
                        double dist = dx*sqrt((i2-i1)*(i2-i1)
                                            + (j2-j1)*(j2-j1)
                                            + (k2-k1)*(k2-k1));
                        //std::cout << i2 << " " << j2 << " " << k2 << " " << std::endl;
                        if(i2<0) i2 += nx;
                        if(j2<0) j2 += ny;
                        if(k2<0) k2 += nz;
                        if(i2>(nx-1)) i2 -= nx;
                        if(j2>(ny-1)) j2 -= ny;
                        if(k2>(nz-1)) k2 -= nz;
                        //arma::vec r2 = {dx*i2, dx*j2, dx*k2};

                        int index = int(dist/delta_x);
                        rhist(index)++;
                        rhist_3d(ind)++;
                        for(int mu=0; mu<3; mu++) 
                        {
                            c_r(index) += xi[t](i1,j1,k1)(mu)*xi[t](i2,j2,k2)(mu);
                            c_r_3d(ind) += xi[t](i1,j1,k1)(mu)*xi[t](i2,j2,k2)(mu);
                        }
                    }
                }
            }
        }
    }
    for(int i=0; i<int(rmax/delta_x)+1; i++) c_r(i) /= (3*rhist(i));
    for(int i=0; i<indices_size; i++) c_r_3d(i) /= (3*rhist_3d(i));

    std::string theName = "spat_corr" + name_add + ".txt";
    std::string the3dName = "spat_corr_3d" + name_add + ".txt";
    print_corr(r, c_r, out_dir, theName);
    print_3d_corr(indices, c_r_3d, dx, out_dir, the3dName);
}

void compute_time_corr(double tmax, double dt, std::vector<arma::field<arma::vec>> &xi, std::string name_add, std::string out_dir)
{
    int nframes = int(xi.size());
    std::cout << "tcorr nframes: " << nframes << std::endl;
    int nx = xi[0].n_rows;
    int ny = xi[0].n_cols;
    int nz = xi[0].n_slices;
    //std::cout << "Dimensions: " << nx << " " << ny << " " << nz << " " << std::endl;
    int nmax = int(tmax/dt);
    arma::vec times(nmax);
    for(int i=0; i<nmax; i++) times(i) = dt*i;
    arma::vec c_t(nmax, arma::fill::zeros);

    std::cout << "Computing time correlation function..." << std::endl;
    for(int t=0; t<nmax; t++)
    {
        if(t%10==0) std::cout << "time difference " << t*dt << std::endl;
        for(int t0=0; t0<nframes-nmax; t0++)
        {
            for(int i=0; i<nx; i++)
            {
                for(int j=0; j<ny; j++)
                {
                    for(int k=0; k<nz; k++)
                    {
                        for(int mu=0; mu<3; mu++)
                        {
                            c_t(t) += xi[t0](i,j,k)(mu)*xi[t0+t](i,j,k)(mu);
                        }                        
                    }
                }
            }
        }
    }
    
    //std::cout << "tcorr normalization: " << 3*nx*ny*nz*(nframes-nmax) << std::endl;
    for(int t=0; t<nmax; t++) c_t(t) /= 3*nx*ny*nz*(nframes-nmax); //3 for 3d

    std::string theName = "time_corr" + name_add + ".txt";
    print_corr(times, c_t, out_dir, theName);
}

void print_corr(arma::vec &coord, arma::vec &corr, std::string out_dir, std::string filename)
{
    assert(coord.size()==corr.size());
    std::ofstream ofile;
    fs::create_directories(out_dir);
    ofile.open(out_dir + "/" + filename);
    for (int i=0; i<int(coord.size()); i++)
    {
        ofile << coord(i) << " " << corr(i) << std::endl;
    }
}

void print_3d_corr(std::vector<std::array<int, 3>> &indices, arma::vec &corr, double dx, std::string out_dir, std::string filename)
{
    std::ofstream ofile;
    fs::create_directories(out_dir);
    ofile.open(out_dir + "/" + filename);
    for(int i=0; i<int(indices.size()); i++)
    {
        ofile << indices[i][0]*dx << " " << indices[i][1]*dx << " " << indices[i][2]*dx << " " << corr(i) << std::endl;
    }
}

arma::vec get_min_image_disp(arma::vec &r1, arma::vec &r2, double Lx, double Ly, double Lz)
{
    arma::vec disp = {0,0,0};
    for (int i=0; i<3; i++) disp[i] = r1[i]-r2[i];

    if (disp[0]<(-0.5*Lx)) disp[0] += Lx;
    if (disp[0]>=(0.5*Lx)) disp[0] -= Lx;

    if (disp[1]<(-0.5*Ly)) disp[1] += Ly;
    if (disp[1]>=(0.5*Ly)) disp[1] -= Ly;

    if (disp[2]<(-0.5*Lz)) disp[2] += Lz;
    if (disp[2]>=(0.5*Lz)) disp[2] -= Lz;

    return disp;
}

double get_min_image_dist(arma::vec &r1, arma::vec &r2, double Lx, double Ly, double Lz)
{
    arma::vec disp = get_min_image_disp(r1, r2, Lx, Ly, Lz);

    double len = 0;
    for (int i=0; i<3; i++) len += disp[i]*disp[i];
    len = sqrt(len);

    return len;
}
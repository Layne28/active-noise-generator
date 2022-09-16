#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <vector>
#include <fstream>
#include <assert.h>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <armadillo>

namespace fs = std::filesystem;
using std::istringstream;

void compute_spatial_corr(double rmax, double delta_x, double Lx, double Ly, double Lz, std::vector<arma::field<arma::vec>> &xi, std::string out_dir);

void compute_time_corr(double tmax, double dt, std::vector<arma::field<arma::vec>> &xi, std::string out_dir);

void print_corr(arma::vec &coord, arma::vec &corr, std::string out_dir, std::string filename);

arma::vec get_min_image_disp(arma::vec &r1, arma::vec &r2, double Lx, double Ly, double Lz);

double get_min_image_dist(arma::vec &r1, arma::vec &r2, double Lx, double Ly, double Lz);

#endif
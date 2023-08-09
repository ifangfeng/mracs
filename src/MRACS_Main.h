#include<fftw3.h>
#include<mkl.h>
#include<mkl_dfti.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstring>
#include<string>
#include<vector>
#include<chrono>
#include<random>
#include<cmath>
#include<omp.h>


#ifndef TWOPI
#define TWOPI M_PI*2
#endif

#ifndef PARAMETERS_NUM
#define PARAMETERS_NUM 10
#endif

#ifndef MRACS_MAIN
#define MRACS_MAIN


extern int Resolution;                    // MRA scale parameter J
extern int BaseType;                      // 0 for B_Spline, 1 for Daubechies phi
extern int phiGenus;                      // Daubechies Wavelet Genus n
extern int SampRate;                      // Wavelet Phi sampling rate (points / 1)
extern int KernelFunc;                    // window function, 0:shell, 1:sphere, 2:Gaussian
extern double Radius;                     // window radius R in Mpc/h
extern double SimBoxL;                    // simulation box length in Mpc/h
extern int Threads;                       // number of threads that used
extern int phiSupport;                    // support of scaling function phi
extern int64_t GridLen;                  // side length of MRA frame, == 2^J
extern int64_t GridVol;                  // number of cubes, == (2^J)^3
extern std::string DIREC;
extern std::string RESOL;
extern std::string RADII;
extern std::string GENUS;
extern std::string DataDirec;
extern std::vector<double> phi;
extern double* PowerPhi;

struct Index
{
    int i;
    int j;
    int k;
};

struct Point
{
    double x;
    double y;
    double z;
};

struct Particle
{
    double x;
    double y;
    double z;
    double weight;
};

struct Halo
{
    double x;
    double y;
    double z;
    double Mass;
    double Concentration;
    double Spin;
};

// double x,y,z,m; int M,C,E,S
struct Hinfo
{
    double x, y, z;
    double mass;
    int Mi{0}, Ci{0}, Ei{0}, Si{0}; // index of mass bin, envir bin, concentration and spin bin;
};

struct Offset
{
    double dx;
    double dy;
    double dz;
    Offset(double x, double y, double z)
    {
        dx = x;
        dy = y;
        dz = z;
    }
};

// defined in mracs.cpp
void welcome();
void read_parameter();
double* sfc_offset(std::vector<Particle>& p, Offset v);
double* sfc_grid_coordinate(std::vector<int64_t>& ps);
double* sfc(std::vector<Particle>& p);
std::vector<double> Daubechies_Phi(const int phiGenus);
std::vector<double> B_Spline(const int n, const int sampleRate);
double* Spectrum1(std::vector<double>& v, double k0, double k1, size_t N_k);
double* Spectrum(std::vector<double>& v, double k0, double k1, size_t N_k);
double* PowerSpectrum(std::vector<double>& v, double k0, double k1, size_t N_k);
double* B_Spline_Dual_Power_Spectrum(double m, double k0, double k1, size_t N_k);
double* PowerPhiFunc(const size_t N);
double* symmetryFold_lean(double* wA);
double* windowArray(const double Radius, const double theta);
double* wft(const double Radius, const double theta);
double* wfc(const double Radius, const double theta);
double* convol3d(double* s, double* w, bool DELETE_S);
double* convol_c2r(fftw_complex* sc, double* w);
double* sc_back(fftw_complex* sc, bool DELETE_SC);
fftw_complex* sfc_r2c(double* s, bool DELETE_S);
fftw_complex* hermitian_product(fftw_complex* sc1, fftw_complex* sc2);
double array_sum(double* w, size_t N);
double inner_product(double* v0, double* v1, size_t N);
double *window_Pk(const double Radius, const double theta);
void force_resoluton_J(int j);
void force_kernel_type(int x);
void force_base_type(int a, int n);
double* prj_grid(const double* s, bool DELETE_S);
double* project_value(const double* s, std::vector<Particle>& p0, bool DELETE_S);

#endif
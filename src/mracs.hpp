#include"mracs_primary.hpp"
#include"kernel.hpp"

// read_in_simdata fuction declared
std::vector<Galaxy> read_in_Millennium_Run_galaxy_catalog(const std::string MilleCata);

//global variables for MRACS interface
int Resolution;                    // MRA scale parameter J
int BaseType;                      // 0 for B_Spline, 1 for Daubechies phi
int phiGenus;                      // Daubechies Wavelet Genus n
int SampRate;                      // Wavelet Phi sampling rate (points / 1)
int KernelFunc;                    // window function, 0:shell, 1:sphere, 2:Gaussian
double Radius;                     // window radius R in Mpc/h
double SimBoxL;                    // simulation box length in Mpc/h
int Threads;                       // number of threads that used

int GridLen;                       // side length of MRA frame, == 2^J
int64_t GridNum;                   // number of cubes, == (2^J)^3
std::string DIREC;
std::string RESOL;
std::string RADII;
std::string GENUS;
std::string MilleCata;
std::vector<double> phi;
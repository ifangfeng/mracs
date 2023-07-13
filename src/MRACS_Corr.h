#include"MRACS_Main.h"

double covar_CombinewithKernel(double* pk_plus, double* WinPk, bool DELET_pk_plus);
double Pk_variance_2dRH(double* Pk, const double Radius, const double theta, int FINE);
double* densityCovarianceArray(fftw_complex* sc1,fftw_complex* sc2);
double* densityVarianceArray(fftw_complex* sc);
double* densityPowerFFT(fftw_complex* sc);
double* densityPowerDWT(fftw_complex* sc);
std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl);
std::vector<double> covar_of_data_vector(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R);
double correlation_coefficients(fftw_complex* sc1, fftw_complex* sc2, double* wpk);
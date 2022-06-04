#include<cmath>

#define TWOPI M_PI*2

double WindowFunction_Shell(double R, double ki, double kj, double kk);
double WindowFunction_Sphere(double R, double ki, double kj, double kk);
double WindowFunction_Gaussian(double R, double ki, double kj, double kk);


//shell window function: sin(kR)/kR
double WindowFunction_Shell(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 1.;
    double Phase = TWOPI * k * R;
    return sin(Phase)/Phase;
};

//sphere window function: 4Pi*[sin(kR)-kRcos(kR)]/k^3
double WindowFunction_Sphere(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 4.*M_PI*pow(R,3)/3.;
    double Phase = TWOPI * k * R;
    return (sin(Phase)-Phase*cos(Phase))/(2*M_PI*M_PI*pow(k,3));
};

//Gaussian window function: e^(-(1/2)(kR)^2)
double WindowFunction_Gaussian(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    double Phase = TWOPI * k * R;
    return pow(1/M_E,Phase*Phase/2);
}
#include"MRACS_Kernel.h"

#ifndef TWOPI
#define TWOPI M_PI*2
#endif



// shell window function: sin(2Pi*kR)/(2Pi*kR)
double WindowFunction_Shell(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 1.;
    double Phase = TWOPI * k * R;
    return sin(Phase)/Phase;
}

// sphere window function: 3*[sin(2Pi*kR)-2Pi*kRcos(2Pi*kR)]/(2Pi*kR)^3
double WindowFunction_Sphere(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 1.;
    double Phase = TWOPI * k * R;
    return 3*(sin(Phase)-Phase*cos(Phase))/(pow(Phase,3));
}

// Gaussian window function: e^(-(1/2)(2Pi*kR)^2)
double WindowFunction_Gaussian(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    double Phase = TWOPI * k * R;
    return pow(1/M_E,Phase*Phase/2);
}

// dual-ring window function: J_0(2Pi*kRsin(a)sin(b))*cos(2Pi*kRcos(a)cos(b))
double WindowFunction_Dual_Ring(double R, double theta, double ki, double kj, double kk)
{
    return (std::cyl_bessel_j(0,TWOPI*sin(theta)*R*sqrt(ki*ki+kj*kj)))*cos(TWOPI * R * cos(theta) * kk);
}

// cylinder window function: sin(2Pi*k_z*h/2)/(Pi*k_z*h/2)*J_1(2Pi*k_r*R)/(2Pi*k_r*R)
double WindowFunction_Cylinder(double R, double h, double ki, double kj, double kk)
{
    double k_r = sqrt(ki * ki + kj * kj);
    return (sin(TWOPI*kk*h/2)/(M_PI*kk*h/2)*std::cyl_bessel_j(1,TWOPI*k_r*R)/(TWOPI*k_r*R));
}
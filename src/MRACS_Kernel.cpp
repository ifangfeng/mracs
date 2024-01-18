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

// thick shell(R1,R2) window function with R1 < R2:
// 3/((2Pi*kR_2)^3-(2Pi*kR_1)^3)*[sin(x)-x*cos(x)]^(2Pi*kR_2)_(2Pi*kR_1)
double WindowFunction_TShell(double R1, double R2, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 1.;
    double Phase1 = TWOPI * k * R1;
    double Phase2 = TWOPI * k * R2;
    return 3*(sin(Phase2)-Phase2*cos(Phase2)-sin(Phase1)+Phase1*cos(Phase1))/(pow(Phase2,3)-pow(Phase1,3));
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

// NFW halo fitting window(), extent == 1 is default, 3 or 10 is an option.
double NFW_window(double r_h, double r_s, int extent, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k==0) return 1.;
    
    double c = r_h / r_s;
    double A = 1. / (log(1+c) - c/(1+c));

    double delta = 0.001;
    int L = c*extent / delta;

    double sum{0};
    for(int i = 1; i < L; ++i){
        sum += sin(TWOPI*k*r_s*i*delta)/(TWOPI*k*r_s*i*delta)*(i*delta)/pow(1+i*delta,2);
    }

    return A*delta*sum;
}

// normalize the kernel to 1 if the 'extent' parameter put into effect, not the default choice
double NFW_window_norm(double r_h, double r_s, int extent, double ki, double kj, double kk)
{
    
    double c = r_h / r_s;
    double A = 1. / log(1+c*extent) - c*extent/(1+c*extent);

    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k==0) return 1.;

    double delta = 0.001;
    int L = c*extent / delta;

    double sum{0};
    for(int i = 1; i < L; ++i){
        sum += sin(TWOPI*k*r_s*i*delta)/(TWOPI*k*r_s*i*delta)*(i*delta)/pow(1+i*delta,2);
    }

    return A*delta*sum;
}
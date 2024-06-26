#include"MRACS_Main.h"

double WindowFunction_Shell(double R, double ki, double kj, double kk);
double WindowFunction_Sphere(double R, double ki, double kj, double kk);
double WindowFunction_Gaussian(double R, double ki, double kj, double kk);
double WindowFunction_TShell(double R1, double R2, double ki, double kj, double kk);
double WindowFunction_Dual_Ring(double R, double theta, double ki, double kj, double kk);
double WindowFunction_Cylinder(double R, double h, double ki, double kj, double kk);
double NFW_window(double r_h, double r_s, int extent, double ki, double kj, double kk);
double NFW_window_norm(double r_h, double r_s, int extent, double ki, double kj, double kk);
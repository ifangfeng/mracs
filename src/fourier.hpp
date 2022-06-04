#include<cmath>

#define TWOPI M_PI*2


void offset(int* d, int n);
void offset2(int* b, int n);
void reorder_lean(double* v, int* d, int n);
void reorder(double* v, int n);
void circle(double* c, int N, int x);
void fft_lean(double* v, double* c, int n);
void FFT(double* v, int n, int Forward);

double* rotate(double* v, int L1, int L2, int L3);
void rotate0(double* vo, double* v, int L1, int L2, int L3);
void FFT3D(double* v, int n1, int n2, int n3, int Forward);
void FFT3D_CUBIC(double* v, int n, int Forward);
void r2rFFT3D_CUBIC(double* v, int n, int Forward);
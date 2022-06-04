#include"fourier.hpp"


//FFT written by feng



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



void offset(int* d, int n)
//reorder element displacement for Real Number
{
    int M, M2{1};
    d[0] = 0;
    for(int m = 0; m < n; ++m)
    {
        M2 <<= 1;
        M = M2 >> 1;
        for(int i = 0; i < M; ++i)
            d[i] = d[i]*2;
        for(int i = M; i < M2; ++i)
            d[i] = d[i-M] + 1;
    }
}

void offset2(int* b, int n)
//reorder element displacement for Complex Number
{
    int N{1<<n};
    int d[N];
    offset(d, n);

    int x;
    for(int i = 0; i < N; ++i)
    {
        x = d[i] << 1;
        b[i+i] = x;
        b[i+i+1]= x + 1;
    }
}

void reorder_lean(double* v, int* d, int n)
//reorder complex number sequence in-place
//d[] comes from offset2
{
    const int N {1<<n};
    const int N2{N<<1};

    double vc[N2];
    for(int i = 0; i < N2; ++i)
        vc[i] = v[i];

    for(int i = 0; i < N2; ++i)
        v[d[i]] = vc[i];
}

void reorder(double* v, int n)
//reorder complex number sequence in-place
{
    const int N {1<<n};
    const int N2{N<<1};

    int d[N2];
    offset2(d, n);
    reorder_lean(v, d, n);
}

void circle(double* c, int N, int x)
//triangle function generated for FFT
{
    const int sign = x ? -1 : 1;
    const int N_2 = N/2;
    //const double TWOPI{2*M_PI};
    const double DeltaPhase = TWOPI / N;

    double Phase{0};

    if(N>2)
        for(int i = 0; i < N_2; i += 2)
        {
            c[i]   = cos(Phase);
            c[i+1] = sin(Phase) * sign;
            c[i+N_2] = c[i+1] * (-sign);
            c[i+N_2+1] = c[i] * sign;
            Phase += DeltaPhase;
        }
    else
    {
        c[0] = 1;
        c[1] = 0;
    }
}

void fft_lean(double* v, double* c, int n)
//1-d complex FFT in place
{
    const int N  {1<<n};
    const int N2 {N<<1};

    int M, M2{1}, M4;
    int step{N2}, jstep;
    double  Real, Image;

    for(int m = 0; m < n; ++m)
    {
        M2 = M2 << 1;
        M  = M2 >> 1;        //M == 2^m
        M4 = M2 << 1;

        step >>= 1;    //pow(2,n-m);
        for(int i = 0; i < N2; i += M4)
            for(int j = 0, jstep = 0; j < M2; j += 2)
            {
                Real  = c[jstep] * v[i + j + M2] - c[jstep + 1] * v[i + j + M2 + 1];
                Image = c[jstep + 1] * v[i + j + M2] + c[jstep] * v[i + j + M2 + 1];
                v[i + j + M2] = -Real + v[i + j];
                v[i + j + M2 + 1] = -Image + v[i + j + 1];
                v[i + j] += Real;
                v[i + j + 1] += Image;
                jstep += step;
            }
    }
}


void FFT(double* v, int n, int Forward)
//1-d complex FFT in place
{
    const int N{1<<n};
    double c[N];

    circle(c, N, Forward);
    reorder(v, n);
    fft_lean(v, c, n);
}

double* rotate(double* v, int L1, int L2, int L3)
//rotate v to a new array v_copy, return the new array as <double> pointer
{
    double* vc = new double[2 * L1 * L2 * L3];

    for(int i = 0; i < L1; ++i)
        for(int j = 0; j < L2; ++j)
            for (int k = 0; k < L3; ++k)
            {
                vc[2 * (i * L2 * L3 + j * L3 + k)] = v[2 * (k * L2 * L3 + i * L3 + j)];
                vc[2 * (i * L2 * L3 + j * L3 + k) + 1] = v[2 * (k * L2 * L3 + i * L3 + j) + 1];
            }
  return vc;
}

void rotate0(double* vo, double* v, int L1, int L2, int L3)
//rotate array v, copy to another array v_original
{
    for(int i = 0; i < L1; ++i)
        for(int j = 0; j < L2; ++j)
            for (int k = 0; k < L3; ++k)
            {
                vo[2 * (i * L2 * L3 + j * L3 + k)] = v[2 * (k * L2 * L3 + i * L3 + j)];
                vo[2 * (i * L2 * L3 + j * L3 + k) + 1] = v[2 * (k * L2 * L3 + i * L3 + j) + 1];
            }
}


void FFT3D(double* v, int n1, int n2, int n3, int Forward)
{
    const int L1  {1<<n1};
    const int L2  {1<<n2};
    const int L3  {1<<n3};
    const int LL1 {L1<<1};
    const int LL2 {L2<<1};
    const int LL3 {L3<<1};
    const int LL  {L1 * L2 * L3};
    const int LLT {LL<<1};

    //====================dimension_3 FFT====================
    int d3[LL3];
    offset2(d3, n3);
    for(int i = 0; i < LLT; i += LL3)
        reorder_lean(&v[i], d3, n3);

    double c3[L3];
    circle(c3, L3, Forward);
    for(int i = 0; i < LLT; i += LL3)
        fft_lean(&v[i], c3, n3);

    //====================dimension_2 FFT====================
    //----before do dimension_2 FFT, rotate v[][][] first----
    double* p = rotate(v, L1, L2, L3);

    int d2[LL2];
    offset2(d2, n2);
    for(int i = 0; i < LLT; i += LL2)
        reorder_lean(&p[i], d2, n2);
    
    double c2[L2];
    circle(c2, L2, Forward);
    for(int i = 0; i < LLT; i += LL2)
        fft_lean(&p[i], c2, n2);

    //====================dimension_1 FFT====================
    //----before do dimension_1 FFT, rotate v[][][] again----
    double* q = rotate(p, L3, L1, L2);
    delete[] p;

    int d1[LL1];
    offset2(d1, n1);
    for(int i = 0; i < LLT; i += LL1)
        reorder_lean(&q[i], d1, n1);
    
    double c1[L1];
    circle(c1, L1, Forward);
    for(int i = 0; i < LLT; i += LL1)
        fft_lean(&q[i], c1, n1);
    
    //============back to original array v[][][]=============
    //---we have done all dimesional FFT, but we need to-----
    //---rotate v[][][] last time to get its orignal order---
    rotate0(v, q, L2, L3, L1);
    delete[] q;
}

void FFT3D_CUBIC(double* v, int n, int Forward)
//same as FFT3D but three equal side length
{
    const int L  {1<<n};
    const int L2 {L<<1};
    const int LL {L * L * L};
    const int LLT {LL<<1};

    int d[L2];
    offset2(d, n);
    double c[L];
    circle(c, L, Forward);
    //====================dimension_3 FFT====================
    for(int i = 0; i < LLT; i += L2)
    {
        reorder_lean(&v[i], d, n);
        fft_lean(&v[i], c, n);
    }
    //====================dimension_2 FFT====================
    //----before do dimension_2 FFT, rotate v[][][] first----
    double* p = rotate(v, L, L, L);
    for(int i = 0; i < LLT; i += L2)
    {
        reorder_lean(&p[i], d, n);
        fft_lean(&p[i], c, n);
    }
    //====================dimension_1 FFT====================
    //----before do dimension_1 FFT, rotate v[][][] again----
    double* q = rotate(p, L, L, L);
    delete[] p;
    for(int i = 0; i < LLT; i += L2)
    {
        reorder_lean(&q[i], d, n);
        fft_lean(&q[i], c, n);
    }
    //============back to original array v[][][]=============
    //---we have done all dimesional FFT, but we need to-----
    //---rotate v[][][] last time to get its orignal order---
    rotate0(v, q, L, L, L);
    delete[] q;
}

//fft 3d-real and store only the real part in place
void r2rFFT3D_CUBIC(double* v, int n, int Forward)
{
    const int L {1<<n};
    auto vc = new double [L*L*L*2]();
    for(int i = 0; i < L*L*L; ++i)
        vc[2*i] = v[i];

    FFT3D_CUBIC(vc, n, Forward);

    for(int i = 0; i < L*L*L; ++i)
        v[i] = vc[2*i];
}










/*
double RFFTPower2(double* v, int j)
{
    const int N{pow(2,j)};
    const int N2{2*N};
    vector <double> Z(N2);

    sin(1/2)*v[1][0]

}

void cmplx_Increase_Two(double* z, int n)
{
    int step{0};
    double wR{0}, wI{0};
    double Real{0}, Image{0};
    for(int i = 0; i < n; ++i)
    {
        step = 2 * i;
        wR = sin(1/n);
        wI = cos(1/n);
        Real      = z[step] * wR - z[step+1] * wI;
        Image     = z[step] * wI + z[step+1] * wR;
        z[step+2] = z[step]   - Real;
        z[step+3] = z[step+1] - Image;
        z[step]   += Real;
        z[step+1] += Image;
    }
}
*/
/*
void compose2(double* v; const int n; int m)
{
    const int N{pow(2,n)};
    const int N2{N*2};
    int M{pow(2,m)};
    int M2{M<<1};
    int M4{M<<2};
    int ij{0};
    double Real, Image;
    double* a = new double[M][2];
    for(int i = 0; i < M; ++i)
    {
        a[i][0] = sin(i/M2);
        a[i][1] = cos(i/M2);
    }
    for(int i = 0; i < N2; i += M4)
        for(int j = 0; j < M2; j += 2)
        {
            ij = i + j + M2;
            Real  = a[j>>1][0] * v[ij] - a[j>>1][1] * v[ij+1];
            Image = a[j>>1][1] * v[ij] + a[j>>1][0] * v[ij+1];
            v[ij] = v[ij-M2] - Real;
            v[ij+1] = v[ij-M2+1] - Image;
            v[ij-M2] += Real;
            v[ij-M2+1] += Image;
        }
    delete [] a;
}

void FFT_Length_2_cmplx(double* v)
{
    v[0] = v[0] + v[2];
    v[1] = v[1] + v[3];
    v[0] = v[0] - v[2];
    v[1] = v[1] - v[3];
}

void FFT_Length_4_cmplx(double* v)
{

}

void re_Sequence_real(double* a, int N)
{
    double temp[N/2];
    for(int i = 0; i < N/2; ++i)
    {
        temp[i] = a[2*i+1];
           a[i] = a[2*i];
    }
    for(i = N/2; i < N; ++i)
        a[i] = temp[i-N/2];
}

void re_Sequence_cmplx(double* a, int n)
{
    double temp[n/2];
    for(int i = 0; i < n/2; ++i)
    {
        temp[i] = a[2*i+1];
        a[i] = a[2*i];
    }
    for(i = n/2; i < n; ++i)
        a[i] = temp[i-n/2];
}
*/
/*
void FFT(double* v, int n, int Forward)
{
    const int N = pow(2,n) ;
    const int N2 {N*2};
    int M2{1};
    int M, M4;
    int ij{0};
    reorder(v, n);
    double Real, Image;
    double TWOPI{2*M_PI};
    double deltaPhase, Phase;
    int sign = Forward ? 1 : -1;
    for(int m = 0; m < n; ++m)
    {
        M2 <<= 1;
        M  = M2>>1;     //M == 2^(m+1)
        M4 = M2<<1;
        double a[N>>1][2];
        deltaPhase = TWOPI / M2;
        Phase = 0;
        for(int i = 0; i < M; ++i)
        {
            a[i][0] = cos(Phase);
            a[i][1] = sin(Phase) * sign;
            Phase += deltaPhase;
        }
        for(int i = 0; i < N2; i += M4)
            for(int j = 0; j < M2; j += 2)
            {
                ij = i + j + M2;
                Real  = a[j>>1][0] * v[ij] - a[j>>1][1] * v[ij+1];
                Image = a[j>>1][1] * v[ij] + a[j>>1][0] * v[ij+1];
                v[ij] = v[ij-M2] - Real;
                v[ij+1] = v[ij-M2+1] - Image;
                v[ij-M2] += Real;
                v[ij-M2+1] += Image;
            }
    }
}
*/
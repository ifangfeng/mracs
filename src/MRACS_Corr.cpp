// =======================================================
// This cpp source file is part of MRACS project
// two point correlation function (Auto and Cross-correlation)
// powerspecturm, cross-corrlation coefficience
// =======================================================
#include"MRACS_Corr.h"


// ************************************************************************************
// return the covariance or variance of density field smoothed by a window w(Radius,theta),
// pk_plus is returned by densityCovarianceArray or densityVarianceArray, WinPk is returned by 
// window_Pk when covarianceArray is the case, covar_CombinewithKernel() returns the covariance.
// ************************************************************************************
double covar_CombinewithKernel(double* pk_plus, double* WinPk, bool DELET_pk_plus)
{
    double sum{0};
    #pragma omp parallel for reduction (+: sum)
    for(int64_t i = 1 - GridLen; i < GridLen; ++i)
        for(int64_t j = 1 - GridLen; j < GridLen; ++j)
            for(int64_t k = 1 - GridLen; k < GridLen; ++k)
            {
                sum += pk_plus[(i+GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (j+GridLen-1) * (2*GridLen-1) + (k+GridLen-1)] *
                WinPk[abs(i) * (GridLen+1) * (GridLen+1) + abs(j) * (GridLen+1) + abs(k)];
            }
    double temp = pk_plus[(GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (GridLen-1) * (2*GridLen-1) + (GridLen-1)] * WinPk[0];

    if(DELET_pk_plus) delete pk_plus;

    return sum/temp-1;
}

// ************************************************************************************
// return the raw covariance i.e. no window smoothing of two density fields,
// pk_plus is returned by densityCovarianceArray or densityVarianceArray.
// ************************************************************************************
double covar_Raw(double* pk_plus, bool DELET_pk_plus)
{
    double sum{0};
    #pragma omp parallel for reduction (+: sum)
    for(int64_t i = 1 - GridLen; i < GridLen; ++i)
        for(int64_t j = 1 - GridLen; j < GridLen; ++j)
            for(int64_t k = 1 - GridLen; k < GridLen; ++k)
            {
                sum += pk_plus[(i+GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (j+GridLen-1) * (2*GridLen-1) + (k+GridLen-1)];
            }
    double temp = pk_plus[(GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (GridLen-1) * (2*GridLen-1) + (GridLen-1)];

    if(DELET_pk_plus) delete pk_plus;

    return sum/temp-1;
}


// ************************************************************************************
// FINE is integral finer parameter and takes value from set {0,1,2,3,4,5}
// notice the Maximun finer parameter is limited to 5 
// ************************************************************************************
double Pk_variance_2dRH(double* Pk, const double Radius, const double theta, int FINE)
{
    if(FINE < 0 || 5 < FINE){
        FINE = 0;
        std::cout << "illegal finer parameter, reset as default FINE=0 " << std::endl;
    }
    const int STEP = 1 << FINE;
    const int64_t LEN {(GridLen/2+1)*STEP};
    double fz[LEN], fr[LEN];
    double* Pk_map = new double[LEN*LEN];

    if(KernelFunc == 3){
        for(size_t i = 0; i < LEN; ++i) fz[i] = cos(TWOPI * i * Radius/SimBoxL * cos(theta)/STEP);
        for(size_t i = 0; i < LEN; ++i) fr[i] = std::cyl_bessel_j(0,TWOPI*i*sin(theta)*Radius/SimBoxL/STEP);
    }
    else if(KernelFunc == 4){
        for(size_t i = 0; i < LEN; ++i) fz[i] = sin(M_PI*i*theta/SimBoxL)/(M_PI*i*theta/SimBoxL);fz[0]=1;
        for(size_t i = 0; i < LEN; ++i) fr[i] = std::cyl_bessel_j(1,TWOPI*Radius/SimBoxL*i)/(M_PI*Radius/SimBoxL*i);fr[0]=1;
    }
    else{
        std::cout << "Pk_variance_2dRH: KernelFunc == 3 or. 4 only" << std::endl;
        return 0;
    }
    for(size_t i = 0; i < LEN; ++i)
        for(size_t j = 0; j < LEN; ++j){
            double k = sqrt(i*i + j*j);
            int kc = k;
            double kf = k-kc;
            int index = k/STEP;
            double xl = Pk[index];
            double xr = Pk[index +1];
            double pk = (xr - xl)*kf + xl;
            Pk_map[i*LEN + j] = pk;
        }
    double var{0};
    for(int i = 0; i < LEN; ++i)
        for(int j = 1-LEN; j < LEN; ++j){
            double ii = i;
            var += fr[i]*fz[abs(j)]*Pk_map[i*LEN+abs(j)]*TWOPI*ii/pow(STEP,3);
        }

    return var;    
}


// ************************************************************************************
// given the sfc array, caluculate its raw covariance array i.e. no window smoothing
// ************************************************************************************
double* densityCovarianceArray(fftw_complex* sc1,fftw_complex* sc2)
{
    auto Pk_array = new double[GridVol];
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]
                + sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
            }
    
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k)
            Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
            Pk_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k];

    auto Pk_plus = new double[(2*GridLen-1)*(2*GridLen-1)*(2*GridLen-1)];

    #pragma omp parallel for
    for(int64_t i = 1 - GridLen; i < GridLen; ++i)
        for(int64_t j = 1 - GridLen; j < GridLen; ++j)
            for(int64_t k = 1 - GridLen; k < GridLen; ++k)
            {
                Pk_plus[(i+GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (j+GridLen-1) * (2*GridLen-1) + (k+GridLen-1)] = 
                Pk_array[((GridLen+i)%GridLen) * GridLen * GridLen + ((GridLen+j)%GridLen) * GridLen + ((GridLen+k)%GridLen)] * 
                PowerPhi[abs(i) * (GridLen+1) * (GridLen+1) + abs(j) * (GridLen+1) + abs(k)];
            }
    
    delete Pk_array;
    
    return Pk_plus;
}


// ************************************************************************************
// given the sfc array, caluculate its raw variance array i.e. no window smoothing
// ************************************************************************************
double* densityVarianceArray(fftw_complex* sc)
{
    return densityCovarianceArray(sc,sc);
}


// ************************************************************************************
// particles first assigned to grid using different window function, then we can take
// advantage of FFT to get its Fourier Coefficiencs and average all orientation
// to have P(k) as function of scalar module k, notice that sfc function do the 
// assignment step exactly
// ************************************************************************************
double* densityPowerFFT(fftw_complex* sc, bool DENOISE)
{
    const double npart = sc[0][0];
    const int64_t PkASize {(GridLen/2 + 1) * GridLen * GridLen};
    double* Pk_array = new double[PkASize];

    #pragma omp parallel for
    for(size_t i = 0; i < PkASize; ++i) Pk_array[i] = pow(sc[i][0], 2) + pow(sc[i][1], 2);

    int64_t klen = GridLen/2 + 1;
    int nk[klen]{0};
    double* Pk = new double[klen]();

    for(size_t i = 0; i < GridLen ; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                int ii = i < GridLen/2 ? i : GridLen - i;
                int jj = j < GridLen/2 ? j : GridLen - j;
                int kM = sqrt(ii * ii + jj * jj + k * k);
                if(kM < klen){
                Pk[kM] += Pk_array[i * GridLen * (GridLen/2+1) + j * (GridLen/2+1) + k];
                nk[kM] += 1;}
            }
    delete[] Pk_array;

    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(npart,2) / pow(SimBoxL,3);
        //Pk[i] -= pow(SimBoxL,3) / npart; // poission shot noise
    }
    Pk[0]=0;

    if(DENOISE){
        for(int i = 1; i < klen; ++i) Pk[i] -= pow(SimBoxL,3) / npart;
    }

    return Pk;
}

//=======================================================================================
// calculate density power spectrum using projected density fileds mathematically
//=======================================================================================
double* densityPowerDWT(fftw_complex* sc, bool DENOISE)
{
    const double npart = sc[0][0];
    double* Pk_array = new double[(GridLen/2 + 1) * GridLen * GridLen];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //#pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                Pk_array[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k] = 
                PowerPhi[(i < GridLen/2 ? i : GridLen - i) * (GridLen + 1) * (GridLen + 1) + (j < GridLen/2 ? j : GridLen - j) * (GridLen + 1) + k] *
                (pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2));
            }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Pk power, T = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    int64_t klen = GridLen/2 + 1;
    int nk[klen]{0};
    double* Pk = new double[klen]();
    
    for(size_t i = 0; i < GridLen ; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                int ii = i < GridLen/2 ? i : GridLen - i;
                int jj = j < GridLen/2 ? j : GridLen - j;
                int kM = sqrt(ii * ii + jj * jj + k * k);
                if(kM < klen)
                {
                Pk[kM] += Pk_array[i * GridLen * (GridLen/2+1) + j * (GridLen/2+1) + k];
                nk[kM] += 1;}
            }
    
    delete[] Pk_array;

    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(npart,2) / pow(SimBoxL,3);
        
    }
    Pk[0]=0;

    if(DENOISE){
        for(int i = 1; i < klen; ++i) Pk[i] -= pow(SimBoxL,3) / npart;
    }
   
    return Pk;
}
//Pk[i] -= pow(SimBoxL,3) / npart; // poission shot noise




//=======================================================================================
// cross-correlation coefficients in Fourier space, using only one realization 
//=======================================================================================
std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl)
{
    auto dm_sc = sfc_r2c(sfc(dm),true);
    auto hl_sc = sfc_r2c(sfc(hl),true);

    double mh[GridLen](0);
    double mm[GridLen](0);
    double hh[GridLen](0);
    
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < (GridLen/2 + 1); ++k){
                auto ii = i < GridLen/2 + 1 ? i : GridLen - i;
                auto jj = j < GridLen/2 + 1 ? j : GridLen - j;
                int ll = sqrt(ii * ii + jj * jj + k * k);
                auto l = i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k;

                mh[ll] += dm_sc[l][0] * hl_sc[l][0] + dm_sc[l][1] * hl_sc[l][1];
                mm[ll] += dm_sc[l][0] * dm_sc[l][0] + dm_sc[l][1] * dm_sc[l][1];
                hh[ll] += hl_sc[l][0] * hl_sc[l][0] + hl_sc[l][1] * hl_sc[l][1];
            }
    std::vector<double> cross(GridLen/2 + 1);
    for(int i = 0; i < GridLen/2 + 1; ++i) cross[i] = mh[i]/sqrt(mm[i] * hh[i]);
    fftw_free(dm_sc);
    fftw_free(hl_sc);

    return cross;
}

std::vector<double> fourier_mode_correlation_1rlz(fftw_complex* dm_sc, fftw_complex* hl_sc)
{
    double mh[GridLen](0);
    double mm[GridLen](0);
    double hh[GridLen](0);
    
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < (GridLen/2 + 1); ++k){
                auto ii = i < GridLen/2 + 1 ? i : GridLen - i;
                auto jj = j < GridLen/2 + 1 ? j : GridLen - j;
                int ll = sqrt(ii * ii + jj * jj + k * k);
                auto l = i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k;

                mh[ll] += dm_sc[l][0] * hl_sc[l][0] + dm_sc[l][1] * hl_sc[l][1];
                mm[ll] += dm_sc[l][0] * dm_sc[l][0] + dm_sc[l][1] * dm_sc[l][1];
                hh[ll] += hl_sc[l][0] * hl_sc[l][0] + hl_sc[l][1] * hl_sc[l][1];
            }
    std::vector<double> cross(GridLen/2 + 1);
    for(int i = 0; i < GridLen/2 + 1; ++i) cross[i] = mh[i]/sqrt(mm[i] * hh[i]);

    return cross;
}
//=======================================================================================
// covariance of vector of particle catalogue smoothed with radius R (kenel type specified in global variable "KernelFunc").
// if the vector in size n, then the returned vector has length (n+1)n/2, cov[i,j]==<vpts[i],vpts[j]>
//=======================================================================================
std::vector<double> covar_of_data_vector(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R)
{
    std::vector<double> cov;
    //auto cov = new std::vector<double>;

    auto wpk = window_Pk(R,0);

    std::vector<fftw_complex*> vec_sc; vec_sc.push_back(sfc_r2c(sfc(dm),true));

    for(auto x : vpts) 
        vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk,true));
    
    return cov;
}

//=======================================================================================
// given two density fields in the form of fourier of sfc coeeficients, and 
// smoothing kernel in the form of |w(k)|^2 as returning by window_pk, 
// return the cross correlation coefficients of these two fileds 
// r = <delta1,delta2> / sqrt(<delta1^2> * <delta2^2>)
//=======================================================================================
double correlation_coefficients(fftw_complex* sc1, fftw_complex* sc2, double* wpk)
{
    auto ab = covar_CombinewithKernel(densityCovarianceArray(sc1,sc2),wpk,true);
    auto aa = covar_CombinewithKernel(densityCovarianceArray(sc1,sc1),wpk,true);
    auto bb = covar_CombinewithKernel(densityCovarianceArray(sc2,sc2),wpk,true);

    double cc = ab/sqrt(aa * bb);

    return cc;
}


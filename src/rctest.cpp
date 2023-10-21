
#include"mracs.h"

fftw_complex* sc_delta(fftw_complex* sc_dm, fftw_complex* sc_hl, double amplitude, bool DELETE_SC_hl)
{
    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    sc[0][0] = sc_dm[0][0];
    sc[0][1] = sc_dm[0][1];
    #pragma omp parallel for
    for(size_t i = 1; i < GridLen * GridLen * (GridLen/2 + 1); ++i){
        sc[i][0] = sc_dm[i][0] - amplitude * sc_hl[i][0];
        sc[i][1] = sc_dm[i][1] - amplitude * sc_hl[i][1];
    }

    if(DELETE_SC_hl) fftw_free(sc_hl);
    return sc;
}
fftw_complex* sc_delta_2(fftw_complex* sc_dm, fftw_complex* sc_hl, double amplitude, bool DELETE_SC_hl)
{
    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i){
        sc[i][0] = sc_dm[i][0]/sc_dm[0][0] - amplitude * sc_hl[i][0]/sc_hl[0][0];
        sc[i][1] = sc_dm[i][1]/sc_dm[0][0] - amplitude * sc_hl[i][1]/sc_hl[0][0];
    }

    if(DELETE_SC_hl) fftw_free(sc_hl);
    return sc;
}
double* densityPowerFFT_2(fftw_complex* sc);

int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto ran = default_random_particle(SimBoxL,hl.size());

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_ran= sfc_r2c(sfc(ran),true);

    std::vector<double*> vec_pk;
    vec_pk.push_back(densityPowerFFT(sc_dm));
    vec_pk.push_back(densityPowerFFT(sc_hl));
    vec_pk.push_back(densityPowerFFT(sc_ran));
    vec_pk.push_back(densityPowerFFT_2(sc_ran));



    for(auto pk : vec_pk){
        for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pk[i] << ", ";std::cout << std::endl;
    }

}

double* densityPowerFFT_2(fftw_complex* sc)
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
                if(kM < GridLen/2 +1){
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

    return Pk;
}


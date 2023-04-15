#include"mracs.h"

int main(){
    read_parameter();
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-5.bin");
    auto s = sfc(p);
    auto sc = sfc_r2c(s);
    auto w = wfc(Radius,0);
    auto c = convol_c2r(sc,w);
    auto var_inner = inner_product(c,c,GridVol)/pow(p.size()/GridLen,2)*GridLen-1;

    std::cout << "BF: " << array_sum(s,GridVol) << "\n";
    std::cout << "AF: " << array_sum(c,GridVol) << "\n";

    std::cout << "var_inner: " << var_inner << std::endl;

    //std::cout << "Fourier:" << std::endl;
    //for(int i = 0; i < 10 ; ++i) std::cout << "[" << sc[i][0] << ", " << sc[i][1] << "], "; std::cout << std::endl;
    //std::cout << "Pk:" << std::endl;
    //for(int i = 0; i < 10 ; ++i) std::cout << sqrt(pow(sc[i][0],2) + pow(sc[i][1],2)) << ", "; std::cout << std::endl;
    auto Pk_af = densityPowerDWT(c);
    double var_pk{0};
    for(size_t i = 0; i < GridVol; ++i) var_pk += Pk_af[i]; --var_pk; std::cout << var_pk << ", and " << Pk_af[0] <<"\n";
    std::cout << "var_pk: " << var_pk << std::endl;
    std::cout << "var_inner: " << var_inner << std::endl;
    std::cout << "var_inner/var_pk: " << var_inner/var_pk << std::endl;

    //std::cout << "=====================a: " << BaseType << " , n: " << phiGenus << "==========================" << std::endl;
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << TWOPI*i/SimBoxL             << ", "; std::cout << '\n' << '\n';
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << Pk[i]            << ", "; std::cout << '\n' << '\n';
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << Pk[i]*i*i*i      << ", "; std::cout << '\n' << '\n';

}
/*
    auto sc = sfc_r2c(s);
    const double npart = array_sum(s,GridVol);
    const uint64_t NyquistL {GridLen/2 + 1};
    double* Pk_array = new double[NyquistL * NyquistL * NyquistL];

    #pragma omp parallel for
    for(size_t i = 0; i < NyquistL; ++i)
        for(size_t j = 0; j < NyquistL; ++j)
            for(size_t k = 0; k < NyquistL; ++k)
            {
                Pk_array[i * NyquistL * NyquistL + j * NyquistL + k] = 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2);
            }
    fftw_free(sc);
    
    #pragma omp parallel for
    for(size_t i = 0; i < NyquistL; ++i)
        for(size_t j = 0; j < NyquistL; ++j)
            for(size_t k = 0; k < NyquistL; ++k)
            {
                Pk_array[i * NyquistL * NyquistL + j * NyquistL + k] *=
                PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k]/pow(npart,2);
            }
    return Pk_array;
*/
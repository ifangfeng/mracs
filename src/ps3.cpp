#include"mracs.h"

int main(){
    read_parameter();
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto s = sfc(p);
    auto sc= sfc_r2c(s);
    //std::cout << "numbers of particle: " << p.size() << std::endl << std::setprecision(10) << "array sum: " << array_sum(s,GridVol) << std::endl
    //<< "Fourier 0 frequence: R= "  << sc[0][0] << " , I= " << sc[0][1] << std::endl;
    auto Pk = densityPowerDWT(sc);
    auto Pk2 = densityPowerFFT(sc);
    //auto Pk3 = densityPowerFFT(s);
    //std::cout << "hh5" << "\n";
    //
//
//
    std::cout << "=====================a: " << BaseType << " , n: " << phiGenus << "==========================" << std::endl;
    for(int i = 0; i < GridLen; ++i) std::cout << TWOPI*i/SimBoxL             << ", "; std::cout << '\n' << "Pk:" << '\n';
    for(int i = 0; i < GridLen; ++i) std::cout << Pk[i]            << ", "; std::cout << '\n' << "Pk2: " << '\n';
    for(int i = 0; i < GridLen; ++i) std::cout << Pk2[i]            << ", "; std::cout << '\n' << "Delta_k: " << '\n';
    for(int i = 0; i < GridLen; ++i) std::cout << Pk[i]*i*i*i      << ", "; std::cout << '\n' << "Delta_k2:" << '\n';
    for(int i = 0; i < GridLen; ++i) std::cout << Pk2[i]*i*i*i      << ", "; std::cout << '\n' << "Delta_k3:" << '\n';
    //for(int i = 0; i < GridLen; ++i) std::cout << Pk3[i]*i*i*i      << ", "; std::cout << '\n' ;
//
    //std::cout << (-1)%GridLen << ", " << (-1-GridLen)%GridLen << ", " << (GridLen-1)%GridLen << ", " << (2*GridLen-1)%GridLen;

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
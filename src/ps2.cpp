#include"mracs.h"

int main(){
    read_parameter();
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto s = sfc(p);
    auto sc = sfc_r2c(s);

    const int NumH{9}, NumD{9}; // shape: H/D; 
    const double Hmin{10}, Dmin{10}, Hmax{100}, Dmax{100};
    std::vector<double> Hvec, Dvec, var_inner, var_pk, var_RE, xvec, yvec;
    for(int i = 0; i <= NumH; ++i) Hvec.push_back(Hmin + i*(Hmax-Hmin)/NumH);
    for(int i = 0; i <= NumD; ++i) Dvec.push_back(Dmin + i*(Dmax-Dmin)/NumD);

std::chrono::steady_clock::time_point begin0 = std::chrono::steady_clock::now();

    for(auto H : Hvec)
        for(auto D : Dvec){
            auto w = wfc(D/2,H);
            auto c = convol_c2r(sc,w);
            var_inner.push_back(inner_product(c,c,GridVol)/pow(p.size(),2)*GridVol-1);
        }

std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    auto pk_plus = densityVarianceArray(sc);
    for(auto H : Hvec)
        for(auto D : Dvec){
            var_pk.push_back(var_CombinewithKernel(pk_plus,D/2,H));
        }

std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    //std::cout << "BF: " << array_sum(s,GridVol) << "\n";
    //std::cout << "AF: " << array_sum(c,GridVol) << "\n";
    //std::cout << "var_inner: " << var_inner << std::endl;

    for(int i = 0; i < var_inner.size(); ++i) var_RE.push_back(var_pk[i]/var_inner[i]-1);
    std::cout << "var_k: " << "\n";for(auto x : var_pk) std::cout << x << ", ";std::cout << std::endl;
    std::cout << "var_x: " << "\n";for(auto x : var_inner) std::cout << x << ", ";std::cout << std::endl;

    std::cout << "MRACS x-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin1 - begin0).count() << "[ms]" << std::endl;
    std::cout << "MRACS k-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin2 - begin1).count() << "[ms]" << std::endl;
    std::cout << "Relative Error (var_k/var_x -1):" << std::endl;
    auto ss = std::cout.precision();
    std::cout << "H vs. D ";for(auto x : Dvec) std::cout << x << ", "; std::cout << std::endl;
    for(int i = 0; i < Hvec.size(); ++i){
        std::cout << std::setprecision(3) << std::defaultfloat << Hvec[i] << ": ";
        for(int j = 0; j < Dvec.size(); ++j) std::cout << std::setprecision(2) << std::scientific << var_RE[i*Dvec.size()+j] << ", ";
        std::cout << std::endl;
        }
    
    //std::cout << "var_pk: " << var_pk << std::endl;
    //std::cout << "var_inner: " << var_inner << std::endl;
    //std::cout << "var_inner/var_pk: " << var_inner/var_pk << std::endl;

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
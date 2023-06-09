#include"mracs.h"

int main(){
    read_parameter();
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl = sfc_r2c(sfc(hl),true);

    int Rmax{20};
    std::vector<double> vec_R, cv_RE, cv_pk, cv_lean; for(int i = 1; i <= Rmax; ++i) vec_R.push_back(i);

     auto pk_plus = densityCovarianceArray(sc_dm,sc_hl);

std::chrono::steady_clock::time_point begin0 = std::chrono::steady_clock::now();
    for(auto r : vec_R){
        auto w = wfc(r,0);
        auto c_dm = convol_c2r(sc_dm,w);
        auto c_hl = convol_c2r(sc_hl,w);
        cv_RE.push_back(GridVol*inner_product(c_dm,c_hl,GridVol)/dm.size()/hl.size() - 1);
    } 
std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

   
    for(auto r : vec_R){
            cv_pk.push_back(var_CombinewithKernel(pk_plus,r,0));
        }

std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();


std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();


    //std::cout << "BF: " << array_sum(s,GridVol) << "\n";
    //std::cout << "AF: " << array_sum(c,GridVol) << "\n";
    //std::cout << "var_inner: " << var_inner << std::endl;

    std::cout << "covar_RE : " << "\n";for(auto x : cv_RE) std::cout << x << ", ";std::cout << std::endl;
    std::cout << "covar_pk : " << "\n";for(auto x : cv_pk) std::cout << x << ", ";std::cout << std::endl;

    //std::cout << "MRACS k-Stat&lean t = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin3 - begin2).count() << "[ms]" << std::endl;
    std::cout << "MRACS x-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin1 - begin0).count() << "[ms]" << std::endl;
    std::cout << "MRACS k-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin2 - begin1).count() << "[ms]" << std::endl;
    
    std::cout << "Relative Error (covar_k/covar_x -1):" << std::endl;
    for(int i = 0; i < cv_pk.size(); ++i) std::cout << cv_pk[i]/cv_RE[i] -1 << ", ";std::cout << std::endl;

    //auto ss = std::cout.precision();
    //std::cout << "H vs. D ";for(auto x : Dvec) std::cout << x << ", "; std::cout << std::endl;
    //for(int i = 0; i < Hvec.size(); ++i){
    //    std::cout << std::setprecision(3) << std::defaultfloat << Hvec[i] << ": ";
    //    for(int j = 0; j < Dvec.size(); ++j) std::cout << std::setprecision(2) << std::scientific << var_RE[i*Dvec.size()+j] << ", ";
    //    std::cout << std::endl;
    //    }
    
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
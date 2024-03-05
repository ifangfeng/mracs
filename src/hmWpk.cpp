#include"mracs.h"

int main(){
    read_parameter();
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto PK = read_in_double("/home/feng/fac/data/Pk_HaloFit_Planck13.bin");

    std::vector<int> vec_R {1,2,3};
    std::ofstream ofs_raw {"output/hmWpk_J"+std::to_string(Resolution)+"raw.txt"};
    //std::ofstream ofs_raw_th {"output/hmWpk_J"+std::to_string(Resolution)+"raw_th.txt"};

    // ----measure---
    force_base_type(1,10);
    auto s = sfc(dm);
    auto sc = sfc_r2c(s,false);

    auto pk = densityPowerFFT(sc,false);
    for(int i = 0; i < GridLen/2+1; ++i) ofs_raw << pk[i] << " "; 
    ofs_raw.close();

    force_kernel_type(2);
    for(auto r : vec_R){
        std::string ofn {"output/hmWpk_J"+std::to_string(Resolution)+"R"+std::to_string(r)+".txt"};
        std::ofstream ofs{ofn};
        auto w = wfc(r,0);
        auto sc_sm = convol_c2c(sc,w);
        auto pk = densityPowerFFT(sc_sm,false);
        for(int i = 0; i < GridLen/2+1; ++i) ofs << pk[i] << " "; 
        delete[] w;
        fftw_free(sc_sm);
    }
    // ----theory----
    //for(auto R : vec_R){
    //    auto win_R = win_theory("sphere",R,0);
    //    std::vector<double> w_tmp(PK.size());
    //    
    //}
    
}
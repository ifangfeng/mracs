#include"mracs.h"

int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);

    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    
    std::vector<double> vec_lambda_opt {4.5, 7.5, 9.25, 16.75, 17.75, 24.0, 23.75};
    std::vector<std::vector<int>> vec_envi;
    for(auto l_th : vec_lambda_opt) {
        auto envi =  web_classify(cxx,hl,l_th);
        vec_envi.push_back(envi);
    }

    delete[] w_gs;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);

    std::vector<double> cc_mul;
    for(size_t i = 0; i < vec_lambda_opt.size(); ++i){
        int Mbin = 1UL<<i;
        std::cout << "Mbin = " << Mbin << "\n";
        double l_th = vec_lambda_opt[i];
        auto mul_vpts  = halo_envi_mass_multi_split(vec_envi[i], hl, Mbin);
        auto sol_mul  = optimal_solution_lean(sc_dm,mul_vpts,wpk,false);
        cc_mul.push_back(sol_mul[0]);        
    }
    std::cout << "cc_mul opt:\n";
    for(auto x : cc_mul) std::cout << x << ", "; std::cout << "\n";
}
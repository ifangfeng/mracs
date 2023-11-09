#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<Particle> hl_n; for(auto x : hl) hl_n.push_back({x.x,x.y,x.z,1.});

    double GSR {1}; // Gaussian smoothing radius
    double THR{15}; // Top-hat smoothing radius

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    
    double lth_Mbin8  = 24.25;
    double lth_Mbin16 = 29.75;

    auto env_8  = web_classify(cxx,hl,lth_Mbin8);
    auto env_16 = web_classify(cxx,hl,lth_Mbin16);

    delete[] w_gs;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);
    
    auto sc_dm = sfc_r2c(sfc(dm),true); std::vector<Particle>().swap(dm);
    auto wpk = window_Pk(THR,0);

    std::vector<double> cc_mul_M, cc_mul_N, cc_M, cc_N;

    auto vpts_8 = halo_envi_mass_multi_split(env_8, hl, 8);
    auto sol_8  = optimal_solution_lean(sc_dm,vpts_8,wpk,false);
    cc_mul_M.push_back(sol_8[0]);
    for(auto x : vpts_8) if(x->size()) std::vector<Particle>().swap(*x);

    auto vpts_16 = halo_envi_mass_multi_split(env_16, hl, 16);
    auto sol_16  = optimal_solution_lean(sc_dm,vpts_16,wpk,false);
    cc_mul_M.push_back(sol_16[0]);
    for(auto x : vpts_16) if(x->size()) std::vector<Particle>().swap(*x);
    
    std::cout << cc_mul_M[0] << ", " << cc_mul_M[1] << "\n";

}
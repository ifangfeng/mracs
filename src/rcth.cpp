#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,4);
    force_kernel_type(2);

    double GSR {1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    
    auto vec_lth = linear_scale_generator(0,50,200,true);
    std::vector<std::vector<int>> vec_env;
    for(auto l_th : vec_lth) {
        auto envi =  web_classify(cxx,hl,l_th);
        vec_env.push_back(envi);
    }

    delete[] w_gs;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);

    double THR{30};
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(THR,0);

    const int LINE{1};
    std::vector<double> cc_mul;
    for(size_t i = 0; i < LINE; ++i){
        int Mbin = 4;//1UL<<i;
        std::cout << "Mbin = " << Mbin << "\n";
        for(auto env : vec_env){
            auto mul_vpts  = halo_envi_mass_multi_split(env, hl, Mbin);
            auto sol_mul  = optimal_solution_lean(sc_dm,mul_vpts,wpk,false);
            cc_mul.push_back(sol_mul[0]);
        }
    }
    for(int i = 0; i < LINE; ++i){
        for(int j = 0; j < vec_lth.size(); ++j)
            std::cout << cc_mul[i * vec_lth.size() + j] << ", ";
        std::cout << std::endl;
    }

}
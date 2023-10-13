#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<double> vec_gsr {8, 5, 3, 2, 1.5, 1, 0.5};
    std::vector<double> vec_lth = linear_scale_generator(0,10,40,true);
    std::vector<std::vector<int>> vec_env;
    std::ofstream ofs {"output/rcgsrTHR15.txt"};

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,4);
    force_kernel_type(2);

    auto sc = sfc_r2c(sfc(dm),true);
    for(auto GSR : vec_gsr){
        auto w_gs = wft(GSR, 0);
        auto cxx = tidal_tensor(sc, w_gs);

        for(auto l_th : vec_lth) {
            auto env =  web_classify(cxx,hl,l_th);
            vec_env.push_back(env);
        }

        delete[] w_gs;
        for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;
    } fftw_free(sc);

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);

    double THR{15};
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(THR,0);

    std::vector<double> cc_env_M, cc_env_N, cc_M;
    // ---N---
    auto cc_N = correlation_coefficients(sc_dm,sfc_r2c(sfc(hl_n),true),wpk);
    // ---env---
    for(size_t i = 0; i < vec_gsr.size(); ++i){
        for(auto env : vec_env){
            auto vpts = halo_envi_match_and_split(env,hl);
            auto sol  = optimal_solution_lean(sc_dm,vpts,wpk,false);
            cc_env_M.push_back(sol[0]);
        }
    }
    for(size_t i = 0; i < vec_gsr.size(); ++i){
        for(auto env : vec_env){
            auto vpts = halo_envi_match_and_split(env,hl_n);
            auto sol  = optimal_solution_lean(sc_dm,vpts,wpk,false);
            cc_env_N.push_back(sol[0]);
        }
    }
    // ---M---
    const int num_Mbin{7};
    for(int i = 0; i < num_Mbin; ++i){
        int Mbin = 1UL<<i; std::cout << "Mbin = " << Mbin << "\n";
        auto vpts = halo_mass_split(hl,Mbin);
        auto sol_M  = optimal_solution_lean(sc_dm,vpts,wpk,false);
        cc_M.push_back(sol_M[0]);
    }
    
    //---output---
    ofs << "----->Mass weighted halo: GSR= "; for(auto x : vec_gsr) ofs << x << ", ";ofs << '\n';
    for(int i = 0; i < vec_gsr.size(); ++i){
        for(int j = 0; j < vec_lth.size(); ++j)
            ofs << cc_env_M[i * vec_lth.size() + j] << ", ";
        ofs << std::endl;
    }
    ofs << "----->Equal weighted halo: GSR= "; for(auto x : vec_gsr) ofs << x << ", ";ofs << '\n';
    for(int i = 0; i < vec_gsr.size(); ++i){
        for(int j = 0; j < vec_lth.size(); ++j)
            ofs << cc_env_N[i * vec_lth.size() + j] << ", ";
        ofs << std::endl;
    }
    ofs << "----->Mass split: Mbin= "; for(int i = 0; i < num_Mbin; ++i) ofs << (1UL<<i) << ", ";ofs << '\n'; 
    for(auto x : cc_M) ofs << x << ", "; ofs << std::endl;
    ofs << "----->default: " << cc_N << std::endl;
    
    

}
#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<Particle> hl_n; for(auto x : hl) hl_n.push_back({x.x,x.y,x.z,1.});

    std::ofstream ofs {"output/rcth_line_fit-THR15.txt"};
    std::string para{"JE="};

    double GSR {1}; // Gaussian smoothing radius
    double THR{15}; // Top-hat smoothing radius

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);
    para+= std::to_string(Resolution);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    
    int lMAX{30};
    auto vec_lth = linear_scale_generator(0,lMAX,120,true);
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
    para+=" , lMAX= " + std::to_string(lMAX) + " , JC="+std::to_string(Resolution);
    para+=" , GSR= "+std::to_string(GSR) + " , THR= "+std::to_string(THR);
    
    auto sc_dm = sfc_r2c(sfc(dm),true); std::vector<Particle>().swap(dm);
    auto wpk = window_Pk(THR,0);

    std::vector<double> cc_mul_M, cc_mul_N, cc_M, cc_N;
    // --- cc_mul_M ----
    const int LINE_M{5};
    for(size_t i = 0; i < LINE_M; ++i){
        int Mbin = 1UL<<i; std::cout << "Mbin = " << Mbin << "\n";
        for(auto env : vec_env){
            auto vpts = halo_envi_mass_multi_split(env, hl, Mbin);
            auto sol  = optimal_solution_lean(sc_dm,vpts,wpk,false);
            cc_mul_M.push_back(sol[0]);
            for(auto x : vpts) if(x->size()) std::vector<Particle>().swap(*x);
        }
    }
    // --- cc_E ---
    //const int LINE_N{3};
    //for(size_t i = 0; i < LINE_N; ++i){
    //    int Nbin = 1UL<<i; std::cout << "Nbin = " << Nbin << "\n";
    //    for(auto env : vec_env){
    //        auto vpts = halo_envi_mass_multi_split(env, hl, Nbin);
    //        for(auto x : vpts) for(auto &y : *x) y.weight = 1.; // make each particle of vpts equal to 1.
    //        auto sol  = optimal_solution_lean(sc_dm,vpts,wpk,false);
    //        cc_mul_N.push_back(sol[0]);
    //        for(auto x : vpts) if(x->size()) std::vector<Particle>().swap(*x);
    //    }
    //}
    // --- cc_M ---
    //const int num_Mbin{3};
    //for(int i = 0; i < num_Mbin; ++i){
    //    int Mbin = 1UL<<i; std::cout << "Mbin = " << Mbin << "\n";
    //    auto vpts = halo_mass_split(hl,Mbin);
    //    auto sol_M  = optimal_solution_lean(sc_dm,vpts,wpk,false);
    //    cc_M.push_back(sol_M[0]);
    //}
    // --- cc_N ---
    //const int num_Nbin{3};
    //for(int i = 0; i < num_Nbin; ++i){
    //    int Nbin = 1UL<<i; std::cout << "Nbin = " << Nbin << "\n";
    //    auto vpts = halo_mass_split(hl,Nbin);
    //    for(auto x : vpts) for(auto &y : *x) y.weight = 1.; // make each particle of vpts equal to 1.
    //    auto sol_N  = optimal_solution_lean(sc_dm,vpts,wpk,false);
    //    cc_N.push_back(sol_N[0]);
    //}

    // ------------------------- output ----------------------------
    ofs << para << '\n';
    ofs << "----->Mass weighted halo:  Mbin= "; for(int i = 0; i < LINE_M; ++i) ofs << (1UL<<i) << ", ";ofs << '\n'; 
    for(int i = 0; i < LINE_M; ++i){
        for(int j = 0; j < vec_lth.size(); ++j)
            ofs << cc_mul_M[i * vec_lth.size() + j] << ", ";
        ofs << std::endl;
    }
    //ofs << "----->Equal weighted halo: Nbin= "; for(int i = 0; i < LINE_N; ++i) ofs << (1UL<<i) << ", ";ofs << '\n'; 
    //for(int i = 0; i < LINE_N; ++i){
    //    for(int j = 0; j < vec_lth.size(); ++j)
    //        ofs << cc_mul_N[i * vec_lth.size() + j] << ", ";
    //    ofs << std::endl;
    //}
    //ofs << "----->Mass split, Mass weighted:  Mbin= "; for(int i = 0; i < num_Mbin; ++i) ofs << (1UL<<i) << ", ";ofs << '\n'; 
    //for(auto x : cc_M) ofs << x << ", "; ofs << std::endl;
    //ofs << "----->Mass split, Equal weighted: Nbin= "; for(int i = 0; i < num_Nbin; ++i) ofs << (1UL<<i) << ", ";ofs << '\n'; 
    //for(auto x : cc_N) ofs << x << ", "; ofs << std::endl;

}
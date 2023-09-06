#include"mracs.h"


int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);
    double GSR {2.1};   // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    
    double l_th{7.5};   // lambda_{th}^{max}=7.5 for Mbin = 2
    auto envi =  web_classify(cxx,hl,l_th);

    delete[] w;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of two mass bin---------
    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);
    int Mbin{2};    // we split halo catalogs with two mass bin each equip four envi parameter

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);

    // ------envi split------
    auto mul_vpts  = halo_envi_mass_multi_split(envi, hl, Mbin);
    auto env_vpts  = halo_envi_match_and_split(envi,hl);
    std::vector<std::vector<Particle>*> vpts_M1Envi, vpts_M2Envi;
    for(int i = 0; auto x : mul_vpts){
        if(i%2==0) vpts_M1Envi.push_back(x);
        else vpts_M2Envi.push_back(x);
        ++i;
    }
    
    auto sol_M1Envi   = optimal_solution_lean(sc_dm,vpts_M1Envi,wpk,false);
    auto sol_M2Envi   = optimal_solution_lean(sc_dm,vpts_M2Envi,wpk,false);
    auto sol_mul  = optimal_solution_lean(sc_dm,mul_vpts,wpk,false);
    auto sol_env  = optimal_solution_lean(sc_dm,env_vpts,wpk,false);

    // ------mass split only-------
    std::cout << "mass-split only:\n";
    auto mass_vpts = halo_mass_split(hl,Mbin);

    std::vector<std::vector<Particle>*> vpts_M1, vpts_M2;
    for(int i = 0; auto x : mass_vpts){
        if(i%2==0) vpts_M1.push_back(x);
        else vpts_M2.push_back(x);
        ++i;
    }
    auto sol_M1   = optimal_solution_lean(sc_dm,vpts_M1,wpk,false);
    auto sol_M2   = optimal_solution_lean(sc_dm,vpts_M2,wpk,false);
    auto sol_mass = optimal_solution_lean(sc_dm,mass_vpts,wpk,false);

    // ----mass bin further split ----
    std::vector<int> v_Mbin{2,4,8};
    std::vector<double> cc_Msplit_M1, cc_Msplit_M2;
    for(auto mbin : v_Mbin){
        auto vpts_M1_Msplit = halo_mass_split(*mass_vpts[0],mbin);
        auto vpts_M2_Msplit = halo_mass_split(*mass_vpts[1],mbin);

        auto sol_M1_Msplit = optimal_solution_lean(sc_dm,vpts_M1_Msplit,wpk,false);
        auto sol_M2_Msplit = optimal_solution_lean(sc_dm,vpts_M2_Msplit,wpk,false);

        cc_Msplit_M1.push_back(sol_M1_Msplit[0]);
        cc_Msplit_M2.push_back(sol_M2_Msplit[0]);
    }

    std::vector<std::vector<Particle>*> vpts_hl{&hl};
    auto sol_non_split = optimal_solution_lean(sc_dm,vpts_hl,wpk,false);

    std::cout << "M1: " << std::setw(8) << sol_M1[0] << ", ----> M1 envi-split: " << std::setw(8) << sol_M1Envi[0] 
    << ", compare with mass split Mbin={2, 4, 8} further: " << cc_Msplit_M1[0] << ", " << cc_Msplit_M1[1] << ", " << cc_Msplit_M1[2] << ", "  << std::endl;
    std::cout << "M2: " << std::setw(8) << sol_M2[0] << ", ----> M2 envi-split: " << std::setw(8) << sol_M2Envi[0] 
    << ", compare with mass split Mbin={2, 4, 8} further: " << cc_Msplit_M2[0] << ", " << cc_Msplit_M2[1] << ", " << cc_Msplit_M2[2] << ", "  << std::endl;
    std::cout << "M: " << std::setw(8) << sol_non_split[0] << std::endl;
    std::cout << "--> M envi-split: " << std::setw(8) << sol_env[0] << std::endl;
    std::cout << "--> M mass-split: " << std::setw(8) << sol_mass[0]<< std::endl;
    std::cout << "--> M E&M multi-split: " << std::setw(8) << sol_mul[0] << std::endl;
}
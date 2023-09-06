// optimal weight and number fraction profile of 4x4 sub-catalogs
#include"mracs.h"


int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environmental classify---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);
    double GSR {2.1};   // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    
    
    // double l_th{9.5};   // lambda_{th}^{max}=7.5 for Mbin = 2
    // auto envi =  web_classify(cxx,hl,l_th);

    std::vector<double> vec_l_th{0, 1, 4.5, 9.5};
    std::vector<std::vector<int>> vec_envi;
    for(auto l_th : vec_l_th) vec_envi.push_back(web_classify(cxx,hl,l_th));

    delete[] w;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of four mass bin---------
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);
    const int Mbin{4}, Ebin{4};    // we split halo catalogs with two mass bin each equip four envi parameter

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);


    double Npart0 = hl.size();
    for(int i = 0; auto envi : vec_envi){
        std::cout << "---->>>> envi_" << i << " with L_th= " << vec_l_th[i] << " :\n";
        auto vpts_env = halo_envi_match_and_split(envi,hl);
        std::cout << "NF:\n";
        for(auto x : vpts_env) std::cout << x->size() / Npart0 << ", "; std::cout << std::endl; 
        auto sol_env  = optimal_solution_lean(sc_dm,vpts_env,wpk,false);
        ++i;
    }
    return 0;
    
    // -------split and solve------
    auto vpts_mul = halo_envi_mass_multi_split(vec_envi[0], hl, Mbin);
    auto sol_mul  = optimal_solution_lean(sc_dm,vpts_mul,wpk,false);

    // -------NF output------
    std::cout << "NF: [M_i x {vd, st, fl, kt}] \n";
    double Npart = hl.size();
    std::vector<double> NF_E(4,0), NF_M(4,0);
    
    for(int i = 0; i < Ebin; ++i){
        for(int j = 0; j < Mbin; ++j){
            NF_E[i] += vpts_mul[i * Mbin + j]->size()/Npart;
        }
    } 
    for(int i = 0; i < Mbin; ++i){
        std::cout << "M" << i << " :";
        for(int j = 0; auto x : vpts_mul) {
            if(j%Mbin==i) {
                std::cout << std::setw(8) << x->size()/Npart << ", ";
                NF_M[i] += x->size()/Npart;
            }
            ++j;
        }
        std::cout << std::setw(8) << NF_M[i] << std::endl;
    }
    for(auto x : NF_E) std::cout << std::setw(8) << x << ", ";
    double NF_E_SUM{0}, NF_M_SUM{0};
    for(auto x : NF_E) NF_E_SUM += x;
    for(auto x : NF_M) NF_M_SUM += x;
    std::cout << "{" << NF_E_SUM << " and " << NF_M_SUM << "}" << std::endl;
    
    // -------NF output------
    std::cout << "W(M,env): [M_i x {vd, st, fl, kt}] \n";
    std::vector<double> w_opt;
    for(int i = 1; i < sol_mul.size(); ++i) w_opt.push_back(sol_mul[i]);
    std::vector<double> W_E(4,0), W_M(4,0);
    for(int i = 0; i < Ebin; ++i){
        for(int j = 0; j < Mbin; ++j){
            W_E[i] += w_opt[i * Mbin + j];
        }
    } 
    for(int i = 0; i < Mbin; ++i){
        std::cout << "M" << i << " :";
        for(int j = 0; auto x : w_opt) {
            if(j%Mbin==i) {
                std::cout << std::setw(8) << x << ", ";
                W_M[i] += x;
            }
            ++j;
        }
        std::cout << std::setw(8) << W_M[i] << std::endl;
    }
    for(auto x : W_E) std::cout << std::setw(8) << x << ", ";
    double W_E_SUM{0}, W_M_SUM{0};
    for(auto x : W_E) W_E_SUM += x;
    for(auto x : W_M) W_M_SUM += x;
    std::cout << "{" << W_E_SUM << " and " << W_M_SUM << "}" << std::endl;

}
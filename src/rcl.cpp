#include"mracs.h"


int main(){
    read_parameter();
    
    //auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");


    // ------cross-correlation of different Lambda_th---------
    std::vector<double> cc_rsol;
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);

    const int Imax{7};

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);

    for(int i = 0; i < Imax; ++i){
        size_t nbin = 1UL << i;
        auto mass_vpts = halo_mass_split(hl,nbin);
        auto sol_mass = optimal_solution_lean(sc_dm,mass_vpts,wpk,false);
        cc_rsol.push_back(sol_mass[0]);
    }
    
    for(auto x : cc_rsol) std::cout << x << ", ";std::cout << "\n";


}
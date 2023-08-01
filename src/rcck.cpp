// density field
#include"mracs.h"


int main(int argc, char** argv){
    read_parameter();
    
    int Mbin{0};
    const int Ebin{4};

    if(argc==1){
        Mbin = 4;
    }
    else {
        int n = std::stoi(argv[1]);
        if ((n >= 1) && (n <= 128))
            Mbin = n;
        else {
            std::cout << "input error, abort\n";
            return 0;
        }
    }

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl_m = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname      {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};

    std::string ofname_hl_n {"output/rcck_hl_n_THR"+RADII+".txt"};
    std::string ofname_hl_m {"output/rcck_hl_m_THR"+RADII+".txt"};
    std::string ofname_rc_M {"output/rcck_M_split_THR" + RADII + ".txt"};
    std::string ofname_rc_ME{"output/rcck_ME_split_THR" + RADII + ".txt"};

    std::ofstream ofshl_n{ofname_hl_n}, ofshl_m{ofname_hl_m},
    ofsrc_M{ofname_rc_M}, ofsrc_ME{ofname_rc_ME};

    std::vector<std::vector<Particle>*> vpts_M, vpts_ME;
    vpts_M = halo_mass_split(hl_m, Mbin);
    vpts_ME = halo_envi_mass_multi_split(ifname, hl_m, Mbin);

    
    // ----reconstruct-------
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl_m = sfc_r2c(sfc(hl_m),true);
    auto sc_hl_n = sfc_r2c(sfc(hl_n),true);

    auto sc_rc_M = optimal_reconstruct(dm,vpts_M,Radius,true);
    auto sc_rc_ME = optimal_reconstruct(dm,vpts_ME,Radius,true);

    auto w = wfc(Radius,0);
    auto s_dm = convol_c2r(sc_dm, w);

    auto ck_n = fourier_mode_correlation_1rlz(sc_dm,sc_hl_n);
    auto ck_m = fourier_mode_correlation_1rlz(sc_dm,sc_hl_m);
    auto ck_M = fourier_mode_correlation_1rlz(sc_dm,sc_rc_M);
    auto ck_ME = fourier_mode_correlation_1rlz(sc_dm,sc_rc_ME);


    for(size_t i = 0; i < ck_n.size(); ++i) ofshl_m << ck_m[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofshl_n << ck_n[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofsrc_M << ck_M[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofsrc_ME << ck_ME[i] << " ";
    
}






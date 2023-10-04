// noise cancelling
#include"mracs.h"


int main(){
    read_parameter();
    int Mbin{4};
    double l_th_opt{9.5};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl_m = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<double> bias{1.16,1.73,1.46,1.47};

    force_resoluton_J(10);
    force_kernel_type(2);
    force_base_type(1,7);

    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    auto env = web_classify(cxx,hl,l_th_opt);

    std::vector<std::vector<Particle>*> vpts_M, vpts_ME;
    vpts_M  = halo_mass_split(hl, Mbin);
    vpts_ME = halo_envi_mass_multi_split(env, hl, Mbin);

    std::vector<double> vec_r{15,30,50,80};

    
    // ----reconstruct-------
    force_resoluton_J(9);
    force_kernel_type(1);
    force_base_type(1,4);

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl_m = sfc_r2c(sfc(hl),true);
    auto sc_hl_n = sfc_r2c(sfc(hl_n),true);

    std::ofstream ofs_hl_n {"output/rcck_Mcut1e13_hl_n.txt"};
    std::ofstream ofs_hl_m {"output/rcck_Mcut1e13_hl_m.txt"};

    auto ck_n = fourier_mode_correlation_1rlz(sc_dm,sc_hl_n);
    auto ck_m = fourier_mode_correlation_1rlz(sc_dm,sc_hl_m);

    for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_m << ck_m[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_n << ck_n[i] << " ";

    for(auto r : vec_r){

        std::ofstream ofs_rc_M {"output/rcck_Mcut1e13_M_split_THRR" + std::to_string(r) + ".txt"};
        std::ofstream ofs_rc_ME{"output/rcck_Mcut1e13_ME_split_THRR" + std::to_string(r) + ".txt"};

        auto sc_rc_M = optimal_reconstruct(dm,vpts_M,r,true);
        auto sc_rc_ME = optimal_reconstruct(dm,vpts_ME,r,true);

        auto ck_M = fourier_mode_correlation_1rlz(sc_dm,sc_rc_M);
        auto ck_ME = fourier_mode_correlation_1rlz(sc_dm,sc_rc_ME);

        for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_M << ck_M[i] << " ";
        for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_ME << ck_ME[i] << " ";

        fftw_free(sc_rc_M);
        fftw_free(sc_rc_ME);
    }

}






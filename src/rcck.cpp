// density field
#include"mracs.h"


int main(){
    read_parameter();

    const int Mbin{4};
    double GSR{1}; // Gaussian smoothing radius
    double lth_opt_ME{17.5}; // optimal lambda_th
    double lth_opt_NE{11.25};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl   = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    auto env_N = web_classify(cxx,hl,lth_opt_NE);
    auto env_M = web_classify(cxx,hl,lth_opt_ME);

    fftw_free(sc);
    delete[] w_gs;
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    std::vector<std::vector<Particle>*> vpts_NE, vpts_ME;
    vpts_NE  = halo_envi_match_and_split(env_N,hl_n);
    vpts_ME = halo_envi_mass_multi_split(env_M, hl, Mbin);

    std::vector<double> vec_r{15,30};

    
    // ----reconstruct-------
    force_resoluton_J(9);
    force_kernel_type(1);
    force_base_type(1,4);

    auto sc_dm = sfc_r2c(sfc(dm),true); 
    auto sc_hl_m = sfc_r2c(sfc(hl),true);
    auto sc_hl_n = sfc_r2c(sfc(hl_n),true);
    std::vector<Particle>().swap(dm);

    std::ofstream ofs_hl_n {"output/rcck_hl_n.txt"};
    std::ofstream ofs_hl_m {"output/rcck_hl_m.txt"};

    auto ck_n = fourier_mode_correlation_1rlz(sc_dm,sc_hl_n);
    auto ck_m = fourier_mode_correlation_1rlz(sc_dm,sc_hl_m);

    for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_m << ck_m[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_n << ck_n[i] << " ";


    std::ofstream ofs_rc_NE{"output/rcck_NE_split.txt"};
    std::ofstream ofs_rc_ME{"output/rcck_ME_split.txt"};

    std::ofstream ofs_rc_NEs{"output/rcck_NE_split_nfw.txt"};
    std::ofstream ofs_rc_MEs{"output/rcck_ME_split_nfw.txt"};

    auto sc_rc_NE = optimal_reconstruct(sc_dm,vpts_NE,r,true);
    auto sc_rc_ME = optimal_reconstruct(sc_dm,vpts_ME,r,true);

    auto ck_NE = fourier_mode_correlation_1rlz(sc_dm,sc_rc_NE);
    auto ck_ME = fourier_mode_correlation_1rlz(sc_dm,sc_rc_ME);

    for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_NE << ck_NE[i] << " ";
    for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_ME << ck_ME[i] << " ";

    fftw_free(sc_rc_NE);
    fftw_free(sc_rc_ME);
    

}






// noise cancelling
#include"mracs.h"

// w is smoothing kernel in fourier space
fftw_complex* sc_delta(fftw_complex* sc_dm, fftw_complex* sc_hl, double amplitude, bool DELETE_SC_hl)
{
    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    sc[0][0] = sc_dm[0][0];
    sc[0][1] = sc_dm[0][1];
    #pragma omp parallel for
    for(size_t i = 1; i < GridLen * GridLen * (GridLen/2 + 1); ++i){
        sc[i][0] = sc_dm[i][0] - amplitude * sc_hl[i][0];
        sc[i][1] = sc_dm[i][1] - amplitude * sc_hl[i][1];
    }

    if(DELETE_SC_hl) fftw_free(sc_hl);
    std::cout << "------delta\n";
    return sc;
}

int main(){
    read_parameter();
    int Mbin{4};
    double l_th_opt{9.5};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl   = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto rd   = default_random_particle(SimBoxL,hl.size());

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
    auto sc_rd = sfc_r2c(sfc(rd),true);

    auto pk_plus_dm = densityVarianceArray(sc_dm);
    auto pk_plus_hl_m = densityVarianceArray(sc_hl_m);
    auto pk_plus_hl_n = densityVarianceArray(sc_hl_n);
    auto pk_plus_rd = densityVarianceArray(sc_rd);

    for(auto r : vec_r){
        std::ofstream ofs_rd {"output/rcnc_rd_split_THRR" + std::to_string(r) + ".txt"};
        std::ofstream ofs_hl_m {"output/rcnc_m_split_THRR" + std::to_string(r) + ".txt"};
        std::ofstream ofs_hl_n {"output/rcnc_n_split_THRR" + std::to_string(r) + ".txt"};
        std::ofstream ofs_rc_M {"output/rcnc_M_split_THRR" + std::to_string(r) + ".txt"};
        std::ofstream ofs_rc_ME{"output/rcnc_ME_split_THRR" + std::to_string(r) + ".txt"};

        auto wpk = window_Pk(r,0);

        auto sc_rc_M = optimal_reconstruct(sc_dm,vpts_M,r,true);
        auto sc_rc_ME = optimal_reconstruct(sc_dm,vpts_ME,r,true);

        double var_dm    = covar_CombinewithKernel(pk_plus_dm,wpk,false);
        double var_rd    = covar_CombinewithKernel(pk_plus_rd,wpk,false);
        double var_hl_m  = covar_CombinewithKernel(pk_plus_hl_m,wpk,false);
        double var_hl_n  = covar_CombinewithKernel(pk_plus_hl_n,wpk,false);
        double var_rc_M  = covar_CombinewithKernel(densityVarianceArray(sc_rc_M),wpk,true);
        double var_rc_ME = covar_CombinewithKernel(densityVarianceArray(sc_rc_ME),wpk,true);

        double b_var_hl_m  = sqrt(var_hl_m/var_dm);
        double b_var_rd    = sqrt(var_rd/var_dm);
        double b_var_hl_n  = sqrt(var_hl_n/var_dm);
        double b_var_rc_M  = sqrt(var_rc_M/var_dm);
        double b_var_rc_ME = sqrt(var_rc_ME/var_dm);

        double cc_hl_m  = correlation_coefficients(sc_hl_m,sc_dm,wpk);
        double cc_rd    = correlation_coefficients(sc_rd,sc_dm,wpk);
        double cc_hl_n  = correlation_coefficients(sc_hl_n,sc_dm,wpk);
        double cc_rc_M  = correlation_coefficients(sc_rc_M,sc_dm,wpk);
        double cc_rc_ME = correlation_coefficients(sc_rc_ME,sc_dm,wpk);

        auto w = wfc(r,0);
        auto pk_hl_m  = densityPowerFFT(sc_delta(sc_dm,sc_hl_m, sc_dm[0][0]/sc_hl_m[0][0] *cc_hl_m/b_var_hl_m,  false));
        auto pk_rd    = densityPowerFFT(sc_delta(sc_dm,sc_rd,   sc_dm[0][0]/sc_rd[0][0]   *cc_rd/b_var_rd,      false));
        auto pk_hl_n  = densityPowerFFT(sc_delta(sc_dm,sc_hl_n, sc_dm[0][0]/sc_hl_n[0][0] *cc_hl_n/b_var_hl_n,  false));
        auto pk_rc_M  = densityPowerFFT(sc_delta(sc_dm,sc_rc_M, sc_dm[0][0]/sc_rc_M[0][0] *cc_rc_M/b_var_rc_M,  false));
        auto pk_rc_ME = densityPowerFFT(sc_delta(sc_dm,sc_rc_ME,sc_dm[0][0]/sc_rc_ME[0][0]*cc_rc_ME/b_var_rc_ME,false));

        delete[] wpk;
        delete[] w;
        fftw_free(sc_rc_M);
        fftw_free(sc_rc_ME);

        for(int i = 0; i < GridLen/2 +1; ++i) ofs_hl_m << pk_hl_m[i] << " ";
        for(int i = 0; i < GridLen/2 +1; ++i) ofs_rd   << pk_rd[i]   << " ";
        for(int i = 0; i < GridLen/2 +1; ++i) ofs_hl_n << pk_hl_n[i] << " ";
        for(int i = 0; i < GridLen/2 +1; ++i) ofs_rc_M << pk_rc_M[i] << " ";
        for(int i = 0; i < GridLen/2 +1; ++i) ofs_rc_ME << pk_rc_ME[i] << " ";

        std::cout << 1./(cc_rd/b_var_rd) << ", " << 1./(cc_hl_n/b_var_hl_n) << ", " << 1./(cc_hl_m/b_var_hl_m) << ", " 
        << 1./(cc_rc_M/b_var_rc_M) << ", " << 1./(cc_rc_ME/b_var_rc_ME) << std::endl;
    }

}






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

    const int Mbin{4};
    double THR{30};
    double GSR{1}; // Gaussian smoothing radius
    double lth_opt_ME{17.5}; // optimal lambda_th
    double lth_opt_NE{11.25};

    std::ofstream ofs_ran {"output/rcnc_ran_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_hl_m {"output/rcnc_m_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_hl_n {"output/rcnc_n_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_rc_NE{"output/rcnc_NE_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_rc_ME{"output/rcnc_ME_split_THRR" + std::to_string(THR) + ".txt"};

    std::vector<std::ofstream*> vec_ofs {&ofs_hl_n, &ofs_hl_m, &ofs_rc_NE, &ofs_rc_ME, &ofs_ran};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl   = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto ran  = default_random_particle(SimBoxL,hl.size());

    std::vector<double> bias{1.16,1.73,1.46,1.47};

    // ------environment sticker---------
    force_resoluton_J(10);
    force_kernel_type(2);
    force_base_type(0,1);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    fftw_free(sc);
    delete[] w_gs;

    auto env_ME = web_classify(cxx,hl,lth_opt_ME);
    auto env_NE = web_classify(cxx,hl,lth_opt_NE);

    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    auto vpts_NE = halo_envi_match_and_split(env_NE,hl_n);
    auto vpts_ME = halo_envi_mass_multi_split(env_ME, hl, Mbin);

    // ----reconstruct-------
    force_resoluton_J(9);
    force_kernel_type(1);
    force_base_type(1,4);

    auto sc_dm   = sfc_r2c(sfc(dm),true);

    std::vector<fftw_complex*> vec_sc;
    vec_sc.push_back(sfc_r2c(sfc(hl_n),true));
    vec_sc.push_back(sfc_r2c(sfc(hl),true));
    vec_sc.push_back(optimal_reconstruct(sc_dm,vpts_NE,THR,true));
    vec_sc.push_back(optimal_reconstruct(sc_dm,vpts_ME,THR,true));

    std::vector<double> vec_a_opt;
    std::vector<double*> vec_pk;

    auto wpk = window_Pk(THR,0);
    for(auto x : vec_sc){
        double a_opt = covar_CombinewithKernel(densityCovarianceArray(x,sc_dm),wpk,true) / 
                        covar_CombinewithKernel(densityVarianceArray(x),wpk,true);
        vec_pk.push_back(densityPowerFFT(sc_delta(sc_dm,x,sc_dm[0][0]/x[0][0] * a_opt, false)));

        vec_a_opt.push_back(a_opt);
    }
    vec_pk.push_back(densityPowerFFT(sfc_r2c(sfc(ran),true)));

    for(int n = 0; n < vec_pk.size(); ++n){
        for(int k = 0; k < GridLen/2 +1; ++k)
            *vec_ofs[n] << vec_pk[n][k] << " ";
        *vec_ofs[n] << std::endl;
    }
    for(auto x : vec_a_opt) std::cout << 1./x << ", "; std::cout << std::endl;

}






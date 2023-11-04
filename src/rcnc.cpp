// noise cancelling
#include"mracs.h"

// w is smoothing kernel in fourier space
fftw_complex* sc_delta(fftw_complex* sc_dm, fftw_complex* sc_hl, double bias, bool DELETE_SC_hl)
{
    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    sc[0][0] = 1;
    sc[0][1] = 0;
    #pragma omp parallel for
    for(size_t i = 1; i < GridLen * GridLen * (GridLen/2 + 1); ++i){
        sc[i][0] = bias * sc_dm[i][0]/sc_dm[0][0] - sc_hl[i][0]/sc_hl[0][0];
        sc[i][1] = bias * sc_dm[i][1]/sc_dm[0][0] - sc_hl[i][1]/sc_hl[0][0];
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

    std::ofstream ofs{"output/rcnc_bias_THR" + std::to_string(THR) + ".txt"};
    
    std::ofstream ofs_hl_n {"output/rcnc_n_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_hl_m {"output/rcnc_m_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_rc_NE{"output/rcnc_NE_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_rc_ME{"output/rcnc_ME_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_ran {"output/rcnc_ran_split_THRR" + std::to_string(THR) + ".txt"};
    std::ofstream ofs_pk_dm{"output/rcnc_pk_dm_split_THRR" + std::to_string(THR) + ".txt"};

    std::vector<std::ofstream*> vec_ofs {&ofs_hl_n, &ofs_hl_m, &ofs_rc_NE, &ofs_rc_ME, &ofs_ran, &ofs_pk_dm};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl   = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto ran  = default_random_particle(SimBoxL,hl.size());

    std::vector<double> bias{1.16,1.73,1.46,1.47};

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);

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
    force_base_type(1,4);
    force_kernel_type(1);

    auto wpk = window_Pk(THR,0);
    auto sc_dm = sfc_r2c(sfc(dm),true);

    trimming_vpts(vpts_NE);
    trimming_vpts(vpts_ME);
    
    std::vector<fftw_complex*> vec_sc_NE, vec_sc_ME;
    vec_sc_NE.push_back(sc_dm); for(auto x : vpts_NE) vec_sc_NE.push_back(sfc_r2c(sfc(*x),true));
    vec_sc_ME.push_back(sc_dm); for(auto x : vpts_ME) vec_sc_ME.push_back(sfc_r2c(sfc(*x),true));

    auto solve_NE = optimal_weight_solver(vec_sc_NE,wpk,true);
    auto solve_ME = optimal_weight_solver(vec_sc_ME,wpk,true);

    reconstruct_with_solve(vpts_NE,solve_NE);
    reconstruct_with_solve(vpts_ME,solve_ME);

    double sum_M1{0}, sum_M2{0};
    double sum_NE1{0}, sum_ME1{0};
    double sum_NE2{0}, sum_ME2{0};
    std::vector<double> vec_expect_noise;
    vec_expect_noise.push_back(pow(SimBoxL,3)/hl.size());

    for(auto x : hl) sum_M1 += x.weight;
    for(auto x : hl) sum_M2 += pow(x.weight,2);
    for(auto x : vpts_NE) for(auto y : *x) sum_NE1 += y.weight;
    for(auto x : vpts_ME) for(auto y : *x) sum_ME1 += y.weight;
    for(auto x : vpts_NE) for(auto y : *x) sum_NE2 += pow(y.weight,2);
    for(auto x : vpts_ME) for(auto y : *x) sum_ME2 += pow(y.weight,2);
    vec_expect_noise.push_back(pow(SimBoxL,3)*sum_M2/pow(sum_M1,2));
    vec_expect_noise.push_back(pow(SimBoxL,3)*sum_NE2/pow(sum_NE1,2));
    vec_expect_noise.push_back(pow(SimBoxL,3)*sum_ME2/pow(sum_ME1,2));

    size_t npts_NE{0}, npts_ME{0};
    for(auto x : vpts_NE) npts_NE += x->size();
    for(auto x : vpts_ME) npts_ME += x->size();

    auto ran_NE = default_random_particle(SimBoxL,npts_NE);
    auto ran_ME = default_random_particle(SimBoxL,npts_ME);

    // weight of vpts_NE to ran_NE
    size_t n{0};
    for(size_t i = 0; i < vpts_NE.size(); ++i){
        for(size_t j = 0; j < vpts_NE[i]->size(); ++j){
            ran_NE[n].weight = (*vpts_NE[i])[j].weight;
            ++n;
        }
    }
    // weight of vpts_ME to ran_ME
    n = 0;
    for(size_t i = 0; i < vpts_ME.size(); ++i){
        for(size_t j = 0; j < vpts_ME[i]->size(); ++j){
            ran_ME[n].weight = (*vpts_ME[i])[j].weight;
            ++n;
        }
    }

    auto sc_ran_NE = sfc_r2c(sfc(ran_NE),true);
    auto sc_ran_ME = sfc_r2c(sfc(ran_ME),true);

    auto pk_ran_NE = densityPowerFFT(sc_ran_NE);
    auto pk_ran_ME = densityPowerFFT(sc_ran_ME);

    
    // ----power spectrum--------
    std::vector<fftw_complex*> vec_sc;
    vec_sc.push_back(sfc_r2c(sfc(hl_n),true));
    vec_sc.push_back(sfc_r2c(sfc(hl),true));
    vec_sc.push_back(reconstruct_with_solve(vec_sc_NE,solve_NE));
    vec_sc.push_back(reconstruct_with_solve(vec_sc_ME,solve_ME));

    std::vector<double> vec_bias;
    std::vector<double*> vec_pk;

    double p_m = covar_CombinewithKernel(densityVarianceArray(sc_dm),wpk,true);
    for(auto x : vec_sc){
        double bias = covar_CombinewithKernel(densityCovarianceArray(x,sc_dm),wpk,true) / p_m;
                        
        vec_pk.push_back(densityPowerFFT(sc_delta(sc_dm, x, bias, false)));

        vec_bias.push_back(bias);
    }
    vec_pk.push_back(densityPowerFFT(sfc_r2c(sfc(ran),true)));
    vec_pk.push_back(densityPowerFFT(sc_dm));

    for(int n = 0; n < vec_pk.size(); ++n){
        for(int k = 0; k < GridLen/2 +1; ++k)
            *vec_ofs[n] << vec_pk[n][k] << " ";
        *vec_ofs[n] << std::endl;
    }
    //--------------------------------

    ofs << "bias: "; for(auto x : vec_bias) ofs << x << ", "; ofs << std::endl;
    ofs << "Noise: "; for(auto x : vec_expect_noise) ofs << x << ", "; ofs << std::endl;
    for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pk_ran_NE[i] << ", "; std::cout << std::endl;
    for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pk_ran_ME[i] << ", "; std::cout << std::endl;
    for(auto x : vec_bias) std::cout << x << ", "; std::cout << std::endl;
    for(auto x :vec_expect_noise) std::cout << x << ", "; std::cout << std::endl;

}




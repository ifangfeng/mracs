// scatter of reconstructed halo vs. dark matter catalogue
#include"mracs.h"


int main(){
    read_parameter();
    
    const int Ebin{4};
    const int Mbin{4};
    double THR{30};
    double GSR{1}; // Gaussian smoothing radius
    double l_th_opt{20}; // optimal lambda_th

    std::string ofname_dm   {"output/rcss_bin"+std::to_string(Mbin)+"_dm_THR"+std::to_string(THR)+".txt"};
    std::string ofname_hl_n {"output/rcss_bin"+std::to_string(Mbin)+"_hl_num_THR"+std::to_string(THR)+".txt"};
    std::string ofname_hl_m {"output/rcss_bin"+std::to_string(Mbin)+"_hl_mass_THR"+std::to_string(THR)+".txt"};
    std::string ofname_rc_M {"output/rcss_bin"+std::to_string(Mbin)+"_M_split_THR" + std::to_string(THR) + ".txt"};
    std::string ofname_rc_ME{"output/rcss_bin"+std::to_string(Mbin)+"_ME_split_THR" + std::to_string(THR) + ".txt"};

    std::ofstream ofsdm{ofname_dm}, ofshl_n{ofname_hl_n}, ofshl_m{ofname_hl_m},
    ofsrc_M{ofname_rc_M}, ofsrc_ME{ofname_rc_ME};

    std::vector<std::vector<Particle>*> vpts_M, vpts_ME;

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_m = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,4);
    force_kernel_type(2);
    
    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    fftw_free(sc);
    
    auto envi =  web_classify(cxx,hl_m,l_th_opt);


    // ----reconstruct-------
    vpts_M = halo_mass_split(hl_m, Mbin);
    vpts_ME = halo_envi_mass_multi_split(envi, hl_m, Mbin);

    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    auto sc_dm   = sfc_r2c(sfc(dm),true);
    auto sc_hl_n = sfc_r2c(sfc(hl_n),true);
    auto sc_hl_m = sfc_r2c(sfc(hl_m),true);
    auto sc_rc_M = optimal_reconstruct(sc_dm,vpts_M,THR,true);
    auto sc_rc_ME = optimal_reconstruct(sc_dm,vpts_ME,THR,true);

    auto w = wfc(THR,0);
    auto s_dm = convol_c2r(sc_dm, w);
    auto s_hl_m = convol_c2r(sc_hl_m, w);
    auto s_hl_n = convol_c2r(sc_hl_n, w);

    auto s_rc_M = convol_c2r(sc_rc_M, w);
    auto s_rc_ME = convol_c2r(sc_rc_ME, w);
    
    // -----sampling------
    double sum_dm{0},sum_hl_m{0},sum_hl_n{0},sum_rc_M{0},sum_rc_ME{0};
    #pragma omp parallel for reduction (+:sum_dm)
    for(size_t i = 0; i < GridVol; ++i) sum_dm += s_dm[i];
    #pragma omp parallel for reduction (+:sum_hl_m)
    for(size_t i = 0; i < GridVol; ++i) sum_hl_m += s_hl_m[i];
    #pragma omp parallel for reduction (+:sum_hl_n)
    for(size_t i = 0; i < GridVol; ++i) sum_hl_n += s_hl_n[i];
    
    #pragma omp parallel for reduction (+:sum_rc_M)
    for(size_t i = 0; i < GridVol; ++i) sum_rc_M += s_rc_M[i];
    #pragma omp parallel for reduction (+:sum_rc_ME)
    for(size_t i = 0; i < GridVol; ++i) sum_rc_ME += s_rc_ME[i];

    size_t N{200};
    auto p0 = generate_random_particle(N,SimBoxL,0);

    auto n_dm = project_value(s_dm,p0,false);
    auto n_hl_m = project_value(s_hl_m,p0,false);
    auto n_hl_n = project_value(s_hl_n,p0,false);
    auto n_rc_M = project_value(s_rc_M,p0,false);
    auto n_rc_ME = project_value(s_rc_ME,p0,false);

    
    for(size_t i = 0; i < N*N*N; ++i) ofsdm << n_dm[i] / sum_dm * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofshl_m << n_hl_m[i] / sum_hl_m * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofshl_n << n_hl_n[i] / sum_hl_n * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofsrc_M << n_rc_M[i] / sum_rc_M * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofsrc_ME << n_rc_ME[i] / sum_rc_ME * GridVol << " ";
    
}






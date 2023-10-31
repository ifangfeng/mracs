// noise cancelling
#include"mracs.h"

int main(){
    read_parameter();

    const int Mbin{4};
    double THR{30};
    double GSR{1}; // Gaussian smoothing radius
    double lth_opt_ME{17.5}; // optimal lambda_th
    double lth_opt_NE{11.25};


    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl   = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

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
    auto pk_NE = densityPowerFFT(sfc_r2c(sfc(ran_NE),true));
    auto pk_ME = densityPowerFFT(sfc_r2c(sfc(ran_ME),true));


    std::vector<double*> vec_pk;

    vec_pk.push_back(pk_NE);
    vec_pk.push_back(pk_ME);


    //--------------------------------
    for(auto x : vec_pk) for(int i = 0; i < GridLen/2+1; ++i) std::cout << x[i] << ", ";std::cout << std::endl;

}




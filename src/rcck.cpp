// density field
#include"mracs.h"


// mass to r_vir
double m2r(double m)
{
    return pow(m,1./3)/49146.1;
}

double m2rs(double m){
    double rs {5.};
    if(m > 1e15) rs = 4.7;
    else if(m > 1e14) rs = 5.4;
    else if(m > 1e13) rs = 6.7;
    else rs = 8.;

    return rs;
}

//=======================================================================================
// window array (convolution kernel)
//=======================================================================================
double* windowArray_NFW(const double r_h, const double r_s)
{
    const double DeltaXi = 1./GridLen;

    auto WindowArray = new double[(GridLen+1) * (GridLen+1) * (GridLen+1)];

    const int refine {10};
    const int arraysize = GridLen * refine * 1.8; // 1.8 is slitly larger than sqrt(3)

    auto fk_nfw = new double[arraysize]; 
    for(int i = 0; i < arraysize; ++i){
        fk_nfw[i] = NFW_window_norm(r_h*GridLen/SimBoxL,r_s*GridLen/SimBoxL,1, 0, 0, i * DeltaXi / refine);
    }

    #pragma omp parallel for
    for(size_t i = 0; i <= GridLen; ++i)
        for(size_t j = 0; j <= GridLen; ++j)
            for(size_t k = 0; k <= GridLen; ++k){
                int index = pow(i*i + j*j + k*k, 0.5) * refine;
                WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = fk_nfw[index];
            }
    return WindowArray;
}

double* wfc_tmp_NFW(const double Radius, const double theta)
{
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto WindowArray = windowArray_NFW(Radius, theta);                          


    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    auto w = symmetryFold_lean(WindowArray);

    return w;
}

double mean_mass(std::vector<Particle>& p){
    if (p.size() == 0) return 0;
    double sum {0};
    for(auto x : p) sum += x.weight;

    return sum/p.size();
}



fftw_complex* optimal_reconstruct_NFW(fftw_complex* sc_dm, std::vector<std::vector<Particle>*>& vpts, double R, bool PRINT)
{
    trimming_vpts(vpts);

    double* wpk = nullptr;

    if(!R)
        wpk = window_Pk(R,0);
    else {
        wpk = new double[(GridLen+1)*(GridLen+1)*(GridLen+1)];
        #pragma omp parallel for
        for(size_t i = 0; i < (GridLen+1)*(GridLen+1)*(GridLen+1); ++i) 
            wpk[i] = 1;
    }

    // -------covariance array-------
    std::vector<fftw_complex*> vec_sc; vec_sc.push_back(sc_dm); 
    for(auto x : vpts) {
        
        auto sc = sfc_r2c(sfc(*x),true);

        auto M = mean_mass(*x);
        auto w = wfc_tmp_NFW(m2r(M),m2rs(M));

        for(size_t i = 0; i < GridLen * GridLen * (GridLen/2+1); ++i) {
            sc[i][0] *= w[i];
            sc[i][1] *= w[i];
        }
        delete[] w;

        vec_sc.push_back(sc);
    }

    // ------solving optimal weight vector------
    auto solve =  optimal_weight_solver(vec_sc,wpk,PRINT);

    // --------reconstruct------------
    return reconstruct_with_solve(vec_sc,solve);
}

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

    //std::ofstream ofs_hl_n {"output/rcck_hl_n_r0.txt"};
    //std::ofstream ofs_hl_m {"output/rcck_hl_m_r0.txt"};
//
    //auto ck_n = fourier_mode_correlation_1rlz(sc_dm,sc_hl_n);
    //auto ck_m = fourier_mode_correlation_1rlz(sc_dm,sc_hl_m);

    //for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_m << ck_m[i] << " ";
    //for(size_t i = 0; i < ck_n.size(); ++i) ofs_hl_n << ck_n[i] << " ";


    //std::ofstream ofs_rc_NE{"output/rcck_NE_split_r0.txt"};
    //std::ofstream ofs_rc_ME{"output/rcck_ME_split_r0.txt"};
//
    //std::ofstream ofs_rc_NEs{"output/rcck_NE_split_nfw_r0.txt"};
    std::ofstream ofs_rc_MEs{"output/rcck_ME_split_nfw_r0.txt"};


    //auto sc_rc_NE = optimal_reconstruct(sc_dm,vpts_NE,Radius,true);
    //auto sc_rc_ME = optimal_reconstruct(sc_dm,vpts_ME,Radius,true);
//
    //auto sc_rc_NEs = optimal_reconstruct_NFW(sc_dm,vpts_NE,Radius,true);
    auto sc_rc_MEs = optimal_reconstruct_NFW(sc_dm,vpts_ME,Radius,true);

    //auto ck_NE = fourier_mode_correlation_1rlz(sc_dm,sc_rc_NE);
    //auto ck_ME = fourier_mode_correlation_1rlz(sc_dm,sc_rc_ME);
    //auto ck_NEs = fourier_mode_correlation_1rlz(sc_dm,sc_rc_NEs);
    auto ck_MEs = fourier_mode_correlation_1rlz(sc_dm,sc_rc_MEs);

    //for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_NE << ck_NE[i] << " ";
    //for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_ME << ck_ME[i] << " ";
//
    //for(size_t i = 0; i < ck_n.size(); ++i) ofs_rc_NEs << ck_NEs[i] << " ";
    for(size_t i = 0; i < ck_MEs.size(); ++i) ofs_rc_MEs << ck_MEs[i] << " ";

    //fftw_free(sc_rc_NE);
    //fftw_free(sc_rc_ME);
    //fftw_free(sc_rc_NEs);
    fftw_free(sc_rc_MEs);
    

}






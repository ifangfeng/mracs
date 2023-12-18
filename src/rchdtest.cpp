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

std::vector<double> mean_mass(std::vector<std::vector<Particle>*>& vpts) {
    std::vector<double> mm;

    for(auto x : vpts) {
        if((*x).size() != 0) {
            double sum = 0; for(auto &y : *x) sum += y.weight;
            mm.push_back(sum/(*x).size());
        }
        else mm.push_back(0);
    }

    return mm;
}

// mass to r_vir
double m2r(double m)
{
    return pow(m,1./3)/49146.1;
}

double* wfc_tmp(const double Radius, const double theta)
{
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto WindowArray = windowArray(Radius, theta);                          

    auto w = symmetryFold_lean(WindowArray);

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    return w;
}

//=======================================================================================
// window array (convolution kernel)
//=======================================================================================
double* windowArray_NFW(const double Radius, const double theta)
{
    const double DeltaXi = 1./GridLen;
    const double RGrid {Radius * GridLen/SimBoxL};

    auto WindowArray = new double[(GridLen+1) * (GridLen+1) * (GridLen+1)];

    #pragma omp parallel for
    for(size_t i = 0; i <= GridLen; ++i)
        for(size_t j = 0; j <= GridLen; ++j)
            for(size_t k = 0; k <= GridLen; ++k)
                WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = NFW_window(RGrid,theta*GridLen/SimBoxL, i * DeltaXi, j * DeltaXi, k * DeltaXi);

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

int main(){
    read_parameter();

    //auto w = wfc_tmp_NFW(2,2/5.7);
    //return 0;

    const int Mbin{4};
    double THR{15};
    double GSR{1}; // Gaussian smoothing radius
    double lth_opt_ME{17.5}; // optimal lambda_th
    double lth_opt_NE{11.25};

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
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

    auto vpts_ME = halo_envi_mass_multi_split(env_ME, hl, Mbin);
    auto vpts_NE = halo_envi_mass_multi_split(env_NE, hl, Mbin);

    trimming_vpts(vpts_NE);
    trimming_vpts(vpts_ME);

    auto mm_ME = mean_mass(vpts_ME);
    auto mm_NE = mean_mass(vpts_NE);

    for(auto x : vpts_NE) for(auto &y : *x) y.weight = 1.; 


    // ----reconstruct-------
    force_resoluton_J(8);
    force_base_type(1,4);
    force_kernel_type(1);

    auto wpk = window_Pk(THR,0);
    auto sc_dm = sfc_r2c(sfc(dm),true);

    std::vector<fftw_complex*> vec_sc_NE, vec_sc_ME, vec_sc_NEs, vec_sc_MEs;
    vec_sc_NE.push_back(sc_dm); for(auto x : vpts_NE) vec_sc_NE.push_back(sfc_r2c(sfc(*x),true));
    vec_sc_ME.push_back(sc_dm); for(auto x : vpts_ME) vec_sc_ME.push_back(sfc_r2c(sfc(*x),true));

    vec_sc_MEs.push_back(sc_dm);
    vec_sc_NEs.push_back(sc_dm);

    for(int i = 0; i < vpts_ME.size(); ++i){
        double r_vir = m2r(mm_ME[i]);
        double c {5.7};
        auto w = wfc_tmp_NFW(r_vir,r_vir/c);
        auto x = convol_c2c(vec_sc_ME[i+1],w);
        vec_sc_MEs.push_back(x);
        delete[] w;
    }
    for(int i = 0; i < vpts_NE.size(); ++i){
        double r_vir = m2r(mm_NE[i]);
        double c {5.7};
        auto w = wfc_tmp_NFW(r_vir,r_vir/c);
        auto x = convol_c2c(vec_sc_NE[i+1],w);
        vec_sc_NEs.push_back(x);
        delete[] w;
    }

    auto solve_NEs = optimal_weight_solver(vec_sc_NEs,wpk,true);
    auto solve_MEs = optimal_weight_solver(vec_sc_MEs,wpk,true);
    auto solve_NE = optimal_weight_solver(vec_sc_NE,wpk,true);
    auto solve_ME = optimal_weight_solver(vec_sc_ME,wpk,true);
    
    std::cout << "r_vir_mean:\n--m: ";
    for(auto x : mm_ME) std::cout << m2r(x) << ", "; std::cout << std::endl << "--n: ";
    for(auto x : mm_NE) std::cout << m2r(x) << ", "; std::cout << std::endl;
    std::cout << "----[Dirac Delta] and [NFW]\n";
    std::cout << "--m:  " << std::setw(10) << solve_ME[0] << ", " << std::setw(10) << solve_MEs[0] << '\n';
    std::cout << "--n:  " << std::setw(10) << solve_NE[0] << ", " << std::setw(10) << solve_NEs[0] << '\n';
}




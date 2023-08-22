#include"mracs.h"


int main(){
    read_parameter();
    
    //auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_kernel_type(2);
    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    
    auto vec_lambda = linear_scale_generator(0,10,40,true);
    std::vector<std::vector<int>> vec_envi;
    for(auto l_th : vec_lambda) {
        auto envi =  web_classify(cxx,hl,l_th);
        vec_envi.push_back(envi);
    }

    delete[] w;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    std::vector<double> cc_mul_rc, mul_data;
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);

    const int nbin{4};
    std::ofstream ofs{"output/rcth4bin.data"};

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);

    for(int i = 0; i < vec_lambda.size(); ++i){
        double l_th = vec_lambda[i];
        std::cout << i << "\n";

        //auto envi_vpts = halo_envi_match_and_split(vec_envi[i], hl);
        auto mul_vpts  = halo_envi_mass_multi_split(vec_envi[i], hl, nbin);
        //auto sol_envi = optimal_solution(dm,envi_vpts,Radius,true);
        auto sol_mul  = optimal_solution_lean(sc_dm,mul_vpts,wpk,false);

        //E_envi_rc.push_back(sqrt(1-pow(sol_envi[0],2)));

        cc_mul_rc.push_back(sol_mul[0]);
        for(auto x : sol_mul) mul_data.push_back(x);
        
    }
    auto mass_vpts = halo_mass_split(hl,nbin);
    auto sol_mass = optimal_solution_lean(sc_dm,mass_vpts,wpk,false);
    std::cout << "mass splited: " << sol_mass[0] << "\n";
    for(auto x : cc_mul_rc) std::cout << x << ", ";std::cout << "\n";

    
    for(auto x : mul_data) ofs << x << " ";
}
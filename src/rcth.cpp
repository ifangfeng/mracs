#include"mracs.h"


int main(){
    read_parameter();
    double GSR {2.1}; // Gaussian smoothing radius
   
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto dm_envi = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    force_kernel_type(2);
    auto sc = sfc_r2c(sfc(dm_envi),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    
    const int nbin{4};

    auto vec_lambda = linear_scale_generator(0,10,40,true);
    std::vector<std::vector<int>> vec_envi;
    for(auto l_th : vec_lambda) {
        auto envi =  web_classify(cxx,hl,l_th);
        vec_envi.push_back(envi);
    }

    std::vector<double> cc_mul_rc;
    force_resoluton_J(7);
    force_kernel_type(1);
    for(int i = 0; i < vec_lambda.size(); ++i){
        double l_th = vec_lambda[i];
        std::cout << i << "\n";

        //auto envi_vpts = halo_envi_match_and_split(vec_envi[i], hl);
        auto mul_vpts  = halo_envi_mass_multi_split(vec_envi[i], hl, nbin);

        //auto sol_envi = optimal_solution(dm,envi_vpts,Radius,true);
        auto sol_mul  = optimal_solution(dm,mul_vpts,Radius,true);

        //E_envi_rc.push_back(sqrt(1-pow(sol_envi[0],2)));
        cc_mul_rc.push_back(sol_mul[0]);
        
    }
    auto mass_vpts = halo_mass_split(hl,nbin);
    auto sol_mass = optimal_solution(dm,mass_vpts,Radius,true);
    std::cout << "mass splited: " << sol_mass[0] << "\n";
    for(auto x : cc_mul_rc) std::cout << x << ", ";std::cout << "\n";


}
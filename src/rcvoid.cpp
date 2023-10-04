#include"mracs.h"


int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);
    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    
    auto vec_lambda = linear_scale_generator(0,25,100,true);
    std::vector<std::vector<int>> vec_envi;
    for(auto l_th : vec_lambda) {
        auto envi =  web_classify(cxx,hl,l_th);
        vec_envi.push_back(envi);
    }

    delete[] w;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(7);
    force_base_type(1,4);
    force_kernel_type(1);

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(Radius,0);


    // ----NF and optimal weight of each envi-type at lth-----
    std::ofstream ofs_w {"output/rcvoid_lth_w.data"};
    std::ofstream ofs_NF {"output/rcvoid_lth_NF.data"};
    std::vector<double> vec_w,vec_NF;
    double Npart = hl.size();
    for(int i = 0; i < vec_lambda.size(); ++i){
        double l_th = vec_lambda[i];
        std::cout << i << "\n";
        auto vpts = halo_envi_match_and_split(vec_envi[i],hl);
        auto sol  = optimal_solution_lean(sc_dm,vpts,wpk,false);

        for(auto x : vpts) vec_NF.push_back(x->size() / Npart);
        for(int i = 1; i < sol.size(); ++i) vec_w.push_back(sol[i]);
    }
    for(auto x : vec_NF) ofs_NF << x << " ";
    for(auto x : vec_w ) ofs_w  << x << " ";

    std::cout << "NF:\n";     for(auto x : vec_NF) std::cout << x << ", ";std::cout << std::endl;
    std::cout << "weight:\n"; for(auto x : vec_w) std::cout << x << ", ";std::cout << std::endl;
}
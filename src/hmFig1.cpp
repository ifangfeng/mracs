#include"mracs.h"
#include<chrono>

int main() {
    read_parameter();
    auto p = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    //auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    
    std::vector<int> vec_resol{8,9,10};
    std::vector<std::ofstream*> vec_ofs;
    std::ofstream ofs_cic{"output/hmFig1_"+RADII+".txt"};
    for(auto x : vec_resol){
        std::ofstream* ofs = new std::ofstream("output/hmFig1_"+RADII+"_J"+std::to_string(x)+".txt");
        vec_ofs.push_back(ofs);
    }

    auto p0 = default_random_particle(SimBoxL,1000);
    for(int i = 0; i < vec_ofs.size(); ++i){
        force_resoluton_J(vec_resol[i]);
        auto s = sfc(p);
        auto w = wfc(Radius,0);
        auto c = convol3d(s,w,true);
        auto n = project_value(c,p0,true);
        const double volume{4./3*M_PI*pow(Radius/SimBoxL*GridLen,3)};
        for(size_t j = 0; j < p0.size(); ++j) *vec_ofs[i] << n[j]*volume << ", ";
        delete[] w;
        delete[] n;
        delete vec_ofs[i];
    }

    auto begin = std::chrono::steady_clock::now();
    auto n = count_in_sphere(Radius,p,p0);
    auto end = std::chrono::steady_clock::now();
    
    for(size_t j = 0; j < p0.size(); ++j) ofs_cic << n[j] << ", ";

    std::cout << "Time difference cic  = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
    << "[ms]" << std::endl;

    delete[] n;
    vec_ofs.clear();
}
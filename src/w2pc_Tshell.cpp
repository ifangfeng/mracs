#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto pk = read_in_double("/home/feng/fac/data/Pk_Zeldovich_Planck13.bin");

    auto vec_r = linear_scale_generator(1,150,40,true);
    std::vector<double> vec_dR{0.1,0.5,1}, xi;

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    force_kernel_type(3);
    for(auto dR : vec_dR) {
        for(auto r : vec_r) {
            auto w = wfc(r,r+dR);
            auto c = convol_c2r(sc,w);
            double dd = inner_product(c,s,GridVol);
            xi.push_back(dd * GridVol/pow(p.size(), 2) - 1);
            delete[] w;
            delete[] c;
        }
    }
    
    for(auto r : vec_r) std::cout << r << ", "; std::cout  << std::endl;

    for(int i = 0; i < vec_dR.size(); ++i) {
        for(int j = 0; j < vec_r.size(); ++j){
            std::cout << xi[i*vec_r.size() + j] << ", "; 
        }
        std::cout << std::endl;
    }
        

}
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");

    auto vec_r = linear_scale_generator(1,150,50,true);
    std::vector<double> xi;

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    force_kernel_type(0);
    for(auto r : vec_r){
        auto w = wfc(r,0);
        auto c = convol_c2r(sc,w);
        double dd = inner_product(c,s,GridVol);
        xi.push_back(dd * GridVol/pow(p.size(), 2) - 1);
        delete[] w;
        delete[] c;
    }
    
    for(auto r : vec_r) 
        std::cout << r << ", "; std::cout  << std::endl;
    for(auto x : xi)
        std::cout << x << ", "; std::cout << std::endl;

}
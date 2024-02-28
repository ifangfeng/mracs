#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");

    auto vec_r = linear_scale_generator(1,1000,50,true);
    std::vector<double> var;

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    force_kernel_type(4);
    for(auto r : vec_r){
        auto w = wfc(r,0);
        auto c = convol_c2r(convol_c2c(sc,w),w);
        var.push_back(inner_product(c,s,GridVol) / GridVol);
        delete[] w;
        delete[] c;
    }

    for(auto r : vec_r) std::cout << r << ", ";std::cout << std::endl;
    std::cout << "var:\n";
    for(auto x : var) 
        std::cout << x << ", "; std::cout  << std::endl;

}
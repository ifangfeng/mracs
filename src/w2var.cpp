#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");

    auto vec_R = linear_scale_generator(1,200,50,true);
    std::vector<double> var;

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    force_kernel_type(4);
    for(auto R : vec_R){
        auto w = wfc(R,0);
        auto c = convol_c2r(sc,w);
        //std::cout << array_sum(c,GridVol) << std::endl;
        var.push_back(inner_product(c,c,GridVol)/pow((p.size()/GridVol),2) / GridVol);
        delete[] w;
        delete[] c;
    }
    
    for(auto R : vec_R) 
        std::cout << R << ", "; std::cout  << std::endl;
    for(auto x : var) std::cout << x << ", "; std::cout << std::endl;

}
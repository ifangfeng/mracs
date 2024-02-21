#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    force_resoluton_J(9);
    

    auto vec_r = linear_scale_generator(1,150,40,true);
    //for(auto x : vec_r) std::cout << x << ", ";std::cout << std::endl;return 0;

    std::vector<double> xi;

    double GSR{15};
    force_kernel_type(2);
    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);
    auto w_gs = wfc(GSR,0);
    auto s_smth = convol3d(convol3d(s,w_gs,false),w_gs,false);

    force_kernel_type(0);
    for(auto r : vec_r){
        auto w = wfc(r,0);
        auto c = convol_c2r(sc,w);
        double dd = inner_product(c,s_smth,GridVol);
        xi.push_back(dd * GridVol/pow(p.size(), 2) - 1);
        delete[] w;
        delete[] c;
    }

    for(auto x : xi) std::cout << x << ", ";std::cout << std::endl;

}
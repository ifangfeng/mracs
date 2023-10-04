#include"mracs.h"

int main(){
    read_parameter();

    auto p0 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    //auto p0 = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");

    force_base_type(1,4);
    force_kernel_type(0);

    auto s0 = sfc(p0);
    auto sc0 = sfc_r2c(s0,false);


    double rmin{1}, rmax{150};
    auto vec_r = linear_scale_generator(rmin,rmax,100,false);

    std::vector<double> xi;
    for(auto i = 0; auto r : vec_r){
        std::cout << "======================R" << i << "=" << r << "\n";
        auto w  = wfc(r,0);

        auto c0 = convol_c2r(sc0,w);

        xi.push_back(inner_product(s0,c0,GridVol)/pow(p0.size(),2)*GridVol -1);

        delete[] w;
        delete[] c0;
        ++i;
    }
    std::cout << "data0:\n"; for(auto x : xi) std::cout << x << ", ";

}
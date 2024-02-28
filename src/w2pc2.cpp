#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");

    auto vec_r = linear_scale_generator(1,150,40,true);
    std::vector<double> vec_R {3,8,15}, xi;

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    auto pk_plus = densityVarianceArray(sc);

    for(auto R : vec_R){
        force_kernel_type(2);
        auto w_smth = wfc(R,0);
        auto sc_smth = convol_c2c(sc,w_smth);
        //auto s_smth = convol3d(convol3d(s,w_smth,false),w_smth,false);

        force_kernel_type(0);
        for(auto r : vec_r){
            
            auto w = wfc(r,0);
            auto c = convol_c2r(sc,w);
            double dd = inner_product(c,s_smth,GridVol);
            xi.push_back(dd * GridVol/pow(p.size(), 2) - 1);
            delete[] w;
            delete[] c;
        }
        delete[] w_smth;
    }
    
    for(auto R : vec_R) 
        std::cout << "smooth R:" << R << ", "; std::cout  << std::endl;
    for(int i = 0; i < vec_R.size(); ++i){
        for(int j = 0; j < vec_r.size(); ++j){
            std::cout << xi[i*vec_r.size()+j] << ", ";
        }std::cout << std::endl;
    }

}
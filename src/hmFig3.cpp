#include"mracs.h"

int main(){
    read_parameter();
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");

    std::ofstream ofs {"output/hmFig3.txt"};
    auto vec_r = linear_scale_generator(1,150,100,false);
    std::vector<double> xi_shell,xi_Tshell;
    auto s  = sfc(p);
    auto sc = sfc_r2c(s,false);
    force_kernel_type(0);
    for(auto r : vec_r){
        auto w = wfc(r,0);
        auto c = convol_c2r(sc,w);
        xi_shell.push_back(inner_product(s,c,GridVol)/pow(p.size(),2)*GridVol -1);
    }
    force_kernel_type(3);
    std::vector<double> vec_thickness {1,3,5};
    for(auto tk : vec_thickness){
        for(auto r : vec_r){
            double R1 = r - tk/2;
            double R2 = r + tk/2;
            if(R1 > 0){
                auto w = wfc(R1,R2);
                auto c = convol_c2r(sc,w);
                xi_Tshell.push_back(inner_product(s,c,GridVol)/pow(p.size(),2)*GridVol -1);
            }
            else{
                auto w = wfc(r-0.25,r+0.25);
                auto c = convol_c2r(sc,w);
                xi_Tshell.push_back(inner_product(s,c,GridVol)/pow(p.size(),2)*GridVol -1);
            }
        }
    }
    for(auto x : vec_r) ofs << x << ", ";
    for(auto x : xi_shell) ofs << x << ", ";
    for(auto x : xi_Tshell) ofs << x << ", ";
    
}
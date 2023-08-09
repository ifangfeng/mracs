#include"mracs.h"


int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-5.bin");
    //auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    auto s = sfc(dm);               std::cout <<std::setprecision(10) << array_sum(s,GridVol) << "\n";
    auto sc = sfc_r2c(s,false);     std::cout <<std::setprecision(10) << array_sum(sc_back(sc,false),GridVol) << "and: " << sc[0][0]<< "\n";

    auto w = wfc(Radius,0);
    auto c = convol_c2r(sc,w);      std::cout <<std::setprecision(10)<< array_sum(c,GridVol) << "\n";

    auto cxx = tidal_tensor(sc,w);
    auto tr = trace_field(cxx);std::cout << "elements: " << tr[0] << ", " << tr[1] << "\n";

    double norm = array_sum(c,GridVol)/GridVol;
    for(size_t i = 0; i < GridVol; ++i){
        c[i] /= norm; c[i] -= 1;
    }
    auto d = prj_grid(c,false);
    double ratio{0};
    size_t count{0};
    for(size_t i = 0; i < GridVol; ++i){
        if(tr[i])
            ratio += d[i]/tr[i];
        else
            count++;
    }std::cout << count << "\n";
    std::cout << ratio/(GridVol-count) << "\n";
}

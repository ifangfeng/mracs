#include"mracs.h"


int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p = default_random_particle(SimBoxL,1000000);
    auto p2= generate_random_particle(100,SimBoxL,0);

    auto sc = sfc_r2c(sfc(p),true);
    auto sc2 = sfc_r2c(sfc(p2),true);
    auto pk = densityPowerFFT(sc);
    auto pk2 = densityPowerFFT(sc2);

    double THR{15};
    auto w = wfc(THR,0);
    auto c = convol3d(sfc(dm),w,false);
    auto sc_c = sfc_r2c(c,false);
    auto pkc= densityPowerFFT(sc_c);


    for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pk[i] << ", ";std::cout << std::endl;
    for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pkc[i] << ", ";std::cout << std::endl;
    for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << pk2[i] << ", ";std::cout << std::endl;
}

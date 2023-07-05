// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"



int main(){
    read_parameter();
    
    //auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto cata = halo_mass_split(hl,12);


}



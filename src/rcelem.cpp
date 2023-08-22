// scatter of reconstructed halo vs. dark matter catalogue
#include"mracs.h"


int main(int argc, char** argv){
    read_parameter();
    
    int Mbin{0};
    const int Ebin{4};

    if(argc==1){
        Mbin = 4;
    }
    else {
        int n = std::stoi(argv[1]);
        if ((n >= 1) && (n <= 128))
            Mbin = n;
        else {
            std::cout << "input error, abort\n";
            return 0;
        }
    }


    std::vector<std::vector<Particle>*> vpts_M, vpts_ME;

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);
    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    fftw_free(sc);
    double l_th_max{7.5}; // optimal lambda_th
    auto envi =  web_classify(cxx,hl,l_th_max);


    // ----reconstruct-------
    vpts_M = halo_mass_split(hl, Mbin);
    vpts_ME = halo_envi_mass_multi_split(envi, hl, Mbin);

    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    auto sc_dm   = sfc_r2c(sfc(dm),true);
    auto sc_hl_m = sfc_r2c(sfc(hl),true);

    auto sc_rc_M = optimal_reconstruct(dm,vpts_M,Radius,true);
    auto sc_rc_ME = optimal_reconstruct(dm,vpts_ME,Radius,true);

    auto w = wfc(Radius,0);
    auto s_dm = convol_c2r(sc_dm, w);
    auto s_hl_m = convol_c2r(sc_hl_m, w);

    auto s_rc_M = convol_c2r(sc_rc_M, w);
    auto s_rc_ME = convol_c2r(sc_rc_ME, w);
    
    
}






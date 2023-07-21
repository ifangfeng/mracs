// covariance of different halo bin 
#include"mracs.h"


int main(int argc, char** argv){
    read_parameter();
    
    int nbin{1};
    if(argc>1){
        int n = std::stoi(argv[1]);
        if(n>=1 && n <= 128){
            nbin = std::stoi(argv[1]);
        }
        else {
            std::cout << "input error, abort\n";
            return 0;
        }
    }

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_uni = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};


    auto envi_vpts = halo_envi_match_and_split(ifname, hl);
    auto mass_vpts = halo_mass_split(hl, nbin);
    auto con_vpts  = halo_envi_mass_concatenate_split(ifname, hl, nbin);
    auto mul_vpts  = halo_envi_mass_multi_split(ifname, hl, nbin);

    // ----reconstruct and check-------
    auto sc_hl_uni = sfc_r2c(sfc(hl_uni),true);
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_dm = sfc_r2c(sfc(dm),true);

    auto sc_envi_rc = optimal_reconstruct(dm,envi_vpts,Radius,true);
    auto sc_mass_rc = optimal_reconstruct(dm,mass_vpts,Radius,true);
    auto sc_con_rc  = optimal_reconstruct(dm,con_vpts,Radius,true);
    auto sc_mul_rc  = optimal_reconstruct(dm,mul_vpts,Radius,true);
    
    auto wpk = window_Pk(Radius,0);

    auto cc_uni = correlation_coefficients(sc_dm,sc_hl_uni,wpk);
    auto cc_hl  = correlation_coefficients(sc_dm,sc_hl,wpk);
    auto cc_envi_rc = correlation_coefficients(sc_dm,sc_envi_rc,wpk);
    auto cc_mass_rc = correlation_coefficients(sc_dm,sc_mass_rc,wpk);
    auto cc_con_rc  = correlation_coefficients(sc_dm,sc_con_rc,wpk);
    auto cc_mul_rc  = correlation_coefficients(sc_dm,sc_mul_rc,wpk);
    
    std::cout << "Cross-correlation coefficient:\n";
    std::cout << "[number picture] default r_n: " << cc_uni      << ", E= " << sqrt(1-pow(cc_uni,2))     << "\n";
    std::cout << "[mass   picture] default r_m: " << cc_hl       << ", E= " << sqrt(1-pow(cc_hl,2))      << "\n";
    std::cout << "---RCST-------envi-------r_m: " << cc_envi_rc  << ", E= " << sqrt(1-pow(cc_envi_rc,2)) << "\n";
    std::cout << "---RCST-------mass-------r_m: " << cc_mass_rc  << ", E= " << sqrt(1-pow(cc_mass_rc,2)) << "\n";
    std::cout << "---RCST----envi + mass---r_m: " << cc_con_rc   << ", E= " << sqrt(1-pow(cc_con_rc,2))  << "\n";
    std::cout << "---RCST----envi * mass---r_m: " << cc_mul_rc   << ", E= " << sqrt(1-pow(cc_mul_rc,2))  << "\n";
    


    
}



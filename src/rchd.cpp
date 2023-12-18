// halo domain
#include"mracs.h"
#include<list>

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-5.bin");
    auto h6 = read_in_Halo_6vector("/data0/MDPL2/halo_Mcut_slice.bin");
    // auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::list<bool> HaloDomain {true,false};
    for(auto HD : HaloDomain){
        std::cout << HD << ", ";
    }std::cout << std::endl;
}

struct cosmology{
    double hubble;
    double Omega_lambda;
    double Omega_m;
    double Omega_b;
    double Sigma_8;
    double n_s;
};


// density sampling of NFW-profile halo catalog, with halo density truncation at TC (times R_vir); 
// DELTA defined as \bar{rho} = M_vir/(4/3*Pi*R_vir^3) = DELTA * rho_crit;
double* NFW_assignment(std::vector<Halo>& halo, double m2r, int TC){
    if (TC < 1 || TC > 5 ) {
        std::cout << "!Warning, halo profile truncation limited to {1,2,3,4,5}, reset TC to 1";
        TC = 1;
    }
    auto s = new double[GridVol](0);
    for(auto x : halo){
        double r_vir = pow(x.Mass,1./3)*m2r;
        int nfw_support = ceil(TC * r_vir / SimBoxL * GridLen);
        
    }
    
}
// Second Order STAtistics (SOSTA)of MRACS code, which include:
// (1) two-point correlation (with or without volume-average) function
// (2) density variance with top-hat filter
#include"mracs.h"

#define R0 0.5           // Mpc/h
#define R1 50.           // Mpc/h
#define NUMTEST 10

int main()
{
    read_parameter();
    std::vector<Particle> p = read_in_DM_3vector(DataDirec);
    // {- uncomment these 3 lines while using Millennium data -}
    // std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    // std::vector<Particle> p; for(Galaxy i : g) p.push_back({i.x, i.y, i.z, 1.});
    // std::vector<Galaxy>().swap(g);
    
    std::vector<double> r_log, xi_r, xi_r_LS, xi_r_ph, xi_r_DP,var_r;

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0;
    for(size_t i = 0; i < p.size(); ++i){
        p0.push_back({u(e), u(e), u(e), p[i].weight});
    }
    
    for(int i = 0; i < NUMTEST; ++i) {
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
    }
    auto s   = sfc(p);
    auto sc  = sfc_r2c(s);
    auto s0  = sfc(p0);
    auto sc0 = sfc_r2c(s0);
    
    // force_kernel_type(0);
    // for(int i = 0; i < NUMTEST; ++i){
    //     auto w = wfc(r_log[i], 0);
    //     auto c = convol_c2r(sc, w);
    //     xi_r.push_back(inner_product(s, c, GridVol) * GridVol/pow(p.size(), 2) - 1);
    //     delete[] c;
    // }

    force_kernel_type(0);
    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        auto c0 = convol_c2r(sc0, w);
        double dd = inner_product(s, c, GridVol);
        double dr = inner_product(s0,c, GridVol);
        double rr = inner_product(s0,c0,GridVol);
        xi_r_DP.push_back(dd/dr-1);
        xi_r_ph.push_back(dd/rr-1);
        xi_r_LS.push_back((dd-2*dr+rr)/rr);
        xi_r.push_back(dd * GridVol/pow(p.size(), 2) - 1);
        delete[] w;
        delete[] c;
    }
    
    force_kernel_type(1);   
    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        var_r.push_back(inner_product(c,c,GridVol)/pow(p.size()*4./3*M_PI*pow(r_log[i]/SimBoxL,3),2)/GridVol-1);
        delete[] c;
    }
    std::cout << "r: " ;
    for(auto x : r_log)  std::cout << x << ", "; std::cout << std::endl << "xi: ";
    for(auto x : xi_r)   std::cout << x << ", "; std::cout << std::endl << "xi_DP: ";
    for(auto x : xi_r_DP)  std::cout << x << ", "; std::cout << std::endl << "xi_ph: ";
    for(auto x : xi_r_ph)  std::cout << x << ", "; std::cout << std::endl << "xi_LS: ";
    for(auto x : xi_r_LS)  std::cout << x << ", "; std::cout << std::endl << "variance: ";
    for(auto x : var_r)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
}

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
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    std::vector<Particle> p, p0; for(Galaxy i : g) p.push_back({i.x, i.y, i.z, i.BulgeMass+i.StellarMass});
    std::vector<double> r_log, xi_r, xi_r_LS, var_r;

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
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
    //     xi_r.push_back(inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1);
    //     delete[] c;
    // }

    force_kernel_type(0);
    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        auto c0 = convol_c2r(sc0, w);
        double dd = inner_product(s, c, GridNum);
        xi_r_LS.push_back((dd-2*inner_product(s,c0,GridNum))/inner_product(s0,c0,GridNum)+1);
        xi_r.push_back(dd * GridNum/pow(p.size(), 2) - 1);
        delete[] w;
        delete[] c;
    }
    
    // force_kernel_type(1);   
    // for( int i = 0; i < NUMTEST; ++i){
    //     auto w = wfc(r_log[i], 0);
    //     auto c = convol_c2r(sc, w);
    //     var_r.push_back(inner_product(c,c,GridNum)/pow(p.size()*4./3*M_PI*pow(r_log[i]/SimBoxL,3),2)/GridNum-1);
    //     delete[] c;
    // }

    for(auto x : r_log) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_r)  std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_r_LS)  std::cout << x << ", "; std::cout << std::endl;
    // for(auto x : var_r) std::cout << sqrt(x) << ", "; std::cout << std::endl;
}
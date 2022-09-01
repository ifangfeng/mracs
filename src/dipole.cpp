#include"mracs.h"

#define RANDOM
#define NUMRAN 9474168
#define NUMTEST 100

int main()
{
    read_parameter();
    #ifdef RANDOM
    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    for(size_t i = 0; i < NUMRAN; ++i) p.push_back({u(e), u(e), u(e), 1});
    #else
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    std::vector<Particle> p; for(auto x : g) p.push_back({x.x, x.y, x.z, x.BulgeMass + x.StellarMass});
    #endif
    std::vector<double> xi_theta, theta;
    for(int i = 1; i < NUMTEST+1; ++i){
        // theta.push_back(acos(1-static_cast<double>(i)/NUMTEST));
        theta.push_back(static_cast<double>(i)/NUMTEST);
    }

    auto s = sfc(p);
    auto sc = sfc_r2c(s);
    force_kernel_type(3);
    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(Radius, theta[i]);
        auto c = convol_c2r(sc, w);
        xi_theta.push_back(inner_product(s, c, GridNum)*GridNum/pow(p.size(), 2) - 1);
        delete[] c;
    }
    // check with shell kernel:
    force_kernel_type(0);
    auto w = wfc(Radius, 0);
    auto c = convol_c2r(sc, w);
    // print out result:
    double sum{0};
    for(auto x : xi_theta) sum += sin(x);
    for(auto x : xi_theta) sum += sin(x);
    for(auto x : theta) std::cout << 2*x/M_PI << ", "; std::cout << std::endl;
    for(auto x : xi_theta) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "ave: " << sum/xi_theta.size() << std::endl;
    std::cout << "shell: " << inner_product(s, c, GridNum)*GridNum/pow(p.size(), 2) - 1 << std::endl;
}

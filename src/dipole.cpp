#include"mracs.h"

#define NUMTEST 100

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p; for(auto x : g) p.push_back({x.x, x.y, x.z, x.BulgeMass + x.StellarMass});
    std::vector<double> xi_theta;

    auto s = sfc(p);
    auto sc = sfc_r2c(s);
    force_kernel_type(3);
    for(int i = 0; i < NUMTEST; ++i)
    {
        auto w = wfc(Radius, acos(static_cast<double>(i)/NUMTEST));//static_cast<double>(i)/NUMTEST*M_PI/2
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
    for(auto x : xi_theta) sum += x;
    for(auto x : xi_theta) std::cout << x << ", ";
    std::cout << std::endl << "ave: " << sum/xi_theta.size() << "\n";
    std::cout << "shell: " << inner_product(s, c, GridNum)*GridNum/pow(p.size(), 2) - 1 << std::endl;
}

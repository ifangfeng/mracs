#include"mracs.hpp"

#define R 10.
#define NUMTEST 1

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(auto x : g) p.push_back({x.x, x.y, x.z, x.BulgeMass + x.StellarMass});

    auto s = sfc(p);
    auto sc = sfc_r2c(s);

    force_kernel_type(3);
    std::vector<double> xi_theta;
    for(int i = 0; i < NUMTEST; ++i)
    {
        auto w = wfc(R, static_cast<double>(i)/NUMTEST*M_PI/2);
        std::cout << w[0] << ", " << w[1] << ", " << w[2] << ", " << w[3] << std::endl;
        auto c = convol_c2r(sc, w);
        xi_theta.push_back(inner_product(s, c, GridNum)*GridNum/pow(p.size(), 2) - 1);
        delete[] c;
    }
    for(auto x : xi_theta) std::cout << x << ", "; std::cout << std::endl;
}

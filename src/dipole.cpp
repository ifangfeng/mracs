#include"mracs.hpp"

#define R 10.
#define TESTPOINTS 20.

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(auto x : g) p.push_back({x.x, x.y, x.z, x.BulgeMass + x.StellarMass});

    auto phi = read_in_phi(phiGenus);
    auto s = scaling_function_coefficients(phi, p);
    auto sc = sfc_r2c(s);

    force_kernel_type(3);
    std::vector<double> xi_theta;
    for(int i = 0; i < TESTPOINTS; ++i)
    {
        auto w = window_function_coefficients(phi, R, static_cast<double>(i)/TESTPOINTS*M_PI/2);
        auto c = convolution_c2r(sc, w);
        xi_theta.push_back(inner_product(s, c, GridNum)*GridNum/pow(p.size(), 2) - 1);
    }
    for(auto x : xi_theta) std::cout << x << ", "; std::cout << std::endl;
}
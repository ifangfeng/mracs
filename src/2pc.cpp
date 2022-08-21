#include"mracs.hpp"

#define R0 10           // Mpc/h
#define R1 50.           // Mpc/h
#define TESTPOINTS 1

int main()
{
    read_parameter();
    force_kernel_type(0);
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy i : g) p.push_back({i.x, i.y, i.z, i.BulgeMass+i.StellarMass});

    auto phi = read_in_phi(phiGenus);
    auto s = scaling_function_coefficients(phi, p);
    auto sc = sfc_r2c(s);

    std::vector<double> r_log, xi_r;
    for(int i = 0; i < TESTPOINTS; ++i)
    {
        std::cout << "|" << std::setw(3) << i << "| th point: ";
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/TESTPOINTS));
        auto w = window_function_coefficients(phi, r_log[i], 0);
        auto c = convolution_c2r(sc, w);
        xi_r.push_back(inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1);
    }
    for(auto x : r_log) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_r)  std::cout << x << ", "; std::cout << std::endl;
}
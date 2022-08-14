#include"mracs.hpp"

#define r0 5
#define r1 50.
#define TESTPOINTS 2

int main()
{
    read_parameter();
    force_kernel_type(1);
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy x : g) p.push_back({x.x, x.y, x.z, x.BulgeMass + x.StellarMass});

    auto phi = read_in_phi(phiGenus);
    auto s = scaling_function_coefficients(phi, p);
    auto sc = sfc_r2c(s);

    std::vector<double> r_log, result;
    for(int i = 0; i < TESTPOINTS; ++i)
    {
        std::cout << "|" << std::setw(3) << i << "| th point: ";
        r_log.push_back(r0 * pow(r1/r0, i/static_cast<double>(TESTPOINTS)));
        auto w = window_function_coefficients(phi, r_log[i], 0);
        auto c = inner_product_c2r(sc, w);
        result.push_back(inner_product(c,c,GridNum)/pow(p.size()*4./3*M_PI*pow(r_log[i]/SimBoxL,3),2)/GridNum-1);
    }
    for(auto x : r_log) std::cout << x << " "; std::cout << std::endl;
    for(auto x : result) std::cout << x << " "; std::cout << std::endl;
}


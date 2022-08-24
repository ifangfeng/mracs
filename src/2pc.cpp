#include"mracs.hpp"

#define R0 0.5           // Mpc/h
#define R1 50.           // Mpc/h
#define NUMTEST 10

int main()
{
    read_parameter();
    force_kernel_type(0);
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy i : g) p.push_back({i.x, i.y, i.z, i.BulgeMass+i.StellarMass});

    auto s = sfc(p);
    auto sc = sfc_r2c(s);

    std::vector<double> r_log, xi_r;
    for(int i = 0; i < NUMTEST; ++i)
    {
        std::cout << "|" << std::setw(3) << i << "| th point: ";
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        xi_r.push_back(inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1);
        delete[] c;
    }
    for(auto x : r_log) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_r)  std::cout << x << ", "; std::cout << std::endl;
}
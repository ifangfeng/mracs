#include"mracs.hpp"

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy i : g)
        p.push_back({i.x, i.y, i.z, 1.});

    auto phi = read_in_phi(phiGenus);
    auto s = scaling_function_coefficients(p, phi, Resolution);
    auto w = window_function_coefficients(phi,Resolution,Radius);
    specialized_convolution_3d(s,w,Resolution);
}
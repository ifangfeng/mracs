// noise cancelling test
// pdf of dm data with different sampling intensities
#include"mracs.h"
#include"kdtree.hpp"

int main()
{
    read_parameter();
    auto p1 = read_in_Halo_4vector("/data0/MDPL2/halo_position.bin");
    std::vector<Particle> p;
    for(auto x : p1) if(x.weight > 2e12) p.push_back({x.x,x.y,x.z,1.});
    std::cout << "number of halo: " << p.size() << std::endl;
    std::vector<Particle>().swap(p1);

    auto p0 = generate_random_particle(1000,SimBoxL,50);

    auto w = wfc(Radius,0);
    auto s = sfc(p);
    auto c = convol3d(s,w);
    double n = 4./3 * M_PI * pow(Radius/SimBoxL,3) * p.size();
    pdf(p0, c, n, 0, 5, 500, "halo_position");
}
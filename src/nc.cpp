// noise cancelling test
// pdf of dm data with different sampling intensities
#include"mracs.h"
#include"kdtree.hpp"

int main()
{
    read_parameter();
    std::vector<std::string> namestr = {"5e-6"};// {"05", "005", "5e-4", "5e-5", "2halo", "5e-6"};
    auto p0 = generate_random_particle(1000,SimBoxL,50);
    auto w = wfc(Radius,0);

    for(auto fname : namestr)
    {
        auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub" + fname + ".bin");
        auto s = sfc(p);
        auto c = convol3d(s,w);
        double n = 4./3 * M_PI * pow(Radius/SimBoxL,3) * p.size();
        pdf(p0, c, n, 0, 5, 20, "dm_sub" + fname);
    }
}


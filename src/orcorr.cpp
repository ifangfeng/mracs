// with-ORiented CORRelation function

#include"mracs.hpp"

#define R 10.           // Mpc/h
#define NUMTEST 20

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy i : g) p.push_back({i.x, i.y, i.z, i.BulgeMass+i.StellarMass});

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0, TWOPI);

    force_kernel_type(0);
    auto s = sfc(p);
    auto w = wfc(R, 0);
    auto c = half_convol(s, w);
    double xi0 = inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1;

    std::vector<Offset> va;
    for(int i = 0; i < NUMTEST; ++i) 
    {
        double a = u(e);
        double b = v(e);
        double c = sqrt(1-a*a);
        va.push_back({c*cos(b)*R, c*sin(b)*R, a*R});
    }

    std::vector<double> xi;
    for(int i = 0; i < NUMTEST; ++i)
    {
        auto o = sfc_offset(p, va[i]);
        xi.push_back(inner_product(s, o, GridNum) * GridNum/pow(p.size(), 2) - 1);
        delete[] o;
    }

    double ave{0}, var{0};
    for(auto x : xi) ave += x;
    for(auto x : xi) var += x*x;
    var = (var-ave*ave/NUMTEST)/(NUMTEST-1);
    ave = ave/NUMTEST;
    for(auto x : xi)  std::cout << x << ", "; std::cout << std::endl;
    std::cout << "ave: " << ave << " ± "<< sqrt(var) << std::endl;
    std::cout << "exp: " << xi0 << std::endl;
}
// with-ORiented CORRelation function

#include"mracs.hpp"

#define R 10.           // Mpc/h
#define TESTPOINTS 20

int main()
{
    read_parameter();
    force_kernel_type(0);
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    std::vector<Particle> p;
    for(Galaxy i : g) p.push_back({i.x, i.y, i.z, i.BulgeMass+i.StellarMass});

    auto phi = read_in_phi(phiGenus);
    auto s = scaling_function_coefficients(phi, p);

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0, TWOPI);

    std::vector<Offset> va;
    for(int i = 0; i < TESTPOINTS; ++i) 
    {
        double a = u(e);
        double b = v(e);
        double c = sqrt(1-a*a);
        va.push_back({c*cos(b)*R, c*sin(b)*R, a*R});
    }

    std::vector<double> xi;
    for(int i = 0; i < TESTPOINTS; ++i)
    {
        auto o = sfc_offset(phi, p, va[i]);
        xi.push_back(inner_product(s, o, GridNum) * GridNum/pow(p.size(), 2) - 1);
        delete[] o;
    }
    
    force_kernel_type(0);
    auto w = window_function_coefficients(phi, R, 0);
    auto c = specialized_convolution_3d(s, w);
    double xi0 = inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1;

    double ave{0}, var{0};
    for(auto x : xi) ave += x;
    for(auto x : xi) var += x*x;
    var = (var-ave*ave/TESTPOINTS)/(TESTPOINTS-1);
    ave = ave/TESTPOINTS;
    for(auto x : xi)  std::cout << x << ", "; std::cout << std::endl;
    std::cout << "ave: " << ave << " ± "<< sqrt(var) << std::endl;
    std::cout << "exp: " << xi0 << std::endl;
}
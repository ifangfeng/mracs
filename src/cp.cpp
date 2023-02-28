#include"mracs.h"

#define NUMTEST  100

int main()
{
    std::ofstream ofsr("output/r.txt");
    std::ofstream ofsprj0("output/prj0.txt");
    std::ofstream ofscic0("output/cic0.txt");
    std::ofstream ofsprj1("output/prj1.txt");
    std::ofstream ofscic1("output/cic1.txt");
    std::ofstream ofsprj2("output/prj2.txt");
    std::ofstream ofscic2("output/cic2.txt");
    if(!ofsr || !ofsprj0 || !ofscic0 || !ofsprj1 || !ofscic1 || !ofsprj2 || !ofscic2 ){
        std::cout << "Abort! opening write out file with error..." << std::endl;
        std::terminate();
    }
    
    read_parameter();
    std::vector<Particle> hl;
    const double M_min {2e12};
    const int NP {3};

    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_position.bin");
    for(size_t i = 0; i < p2.size(); ++i)
        if(p2[i].weight > M_min) 
            hl.push_back({p2[i].x, p2[i].y, p2[i].z, 1.});
    std::vector<Particle>().swap(p2);
    std::cout << "hl: " << hl.size() << std::endl;

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; p0.push_back({127.1, 237.6, 328.6, 1.});p0.push_back({547.1, 337.6, 328.6, 1.});p0.push_back({527.8, 237.3, 398.6, 1.});

    std::vector<double> rhopj0(NUMTEST), rhopj1(NUMTEST), rhopj2(NUMTEST), rhoc0(NUMTEST), rhoc1(NUMTEST), rhoc2(NUMTEST);
    
    const double r0 {5};
    const double r1 {50};
    const double dr {(r1 - r0) / NUMTEST};
    std::vector<double> r(NUMTEST); for(int i = 0; i < NUMTEST; ++i) r[i] = r0 + i * dr;

    force_kernel_type(1);
    auto s  = sfc(hl);
    auto sc = sfc_r2c(s);

    for(int i = 0; i < NUMTEST; ++i)
    {
        auto w = wfc(r[i],0);
        auto c = convol_c2r(sc,w);
        auto temp = project_value(c,p0);
        double volume = 4./3 * M_PI * pow(r[i],3);
        rhopj0[i] = temp[0] / volume;
        rhopj1[i] = temp[1] / volume;
        rhopj2[i] = temp[2] / volume;
        delete[] w;
        delete[] c;
    }

    for(int i = 0; i < NUMTEST; ++i) 
    {
        auto temp = count_in_sphere(r[i],hl,p0);
        double volume = 4./3 * M_PI * pow(r[i],3);
        rhoc0[i] = temp[0] / volume;
        rhoc1[i] = temp[1] / volume;
        rhoc2[i] = temp[2] / volume;
    }

    const double rhobar = hl.size() / pow(SimBoxL,3);
    for(auto i : r) ofsr << i << ", ";
    for(auto i : rhopj0) ofsprj0 << i / rhobar << ", ";
    for(auto i : rhopj1) ofsprj1 << i / rhobar << ", ";
    for(auto i : rhopj2) ofsprj2 << i / rhobar << ", ";
    for(auto i : rhoc0)  ofscic0 << i / rhobar << ", ";
    for(auto i : rhoc1)  ofscic1 << i / rhobar << ", ";
    for(auto i : rhoc2)  ofscic2 << i / rhobar << ", ";


    
}

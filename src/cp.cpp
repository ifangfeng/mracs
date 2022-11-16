#include"mracs.h"

#define NUMRAN 1000

int main()
{
    std::string prjfname {"prj.txt"};
    std::string cicfname {"cic.txt"};
    std::ofstream ofsprj(prjfname);
    std::ofstream ofscic(cicfname);
    if(!ofsprj || !ofscic)
    {
        std::cout << "Abort! opening write out file with error..." << std::endl;
        std::terminate();
    }
    
    read_parameter();
    auto p = read_in_Halo_4vector("/data0/BigMDPL/BigMDPL_halo.bin");

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; for(size_t i = 0; i < NUMRAN; ++i) p0.push_back({u(e), u(e), u(e), 1.});

    force_kernel_type(1);
    auto s  = sfc(p);
    auto sc = sfc_r2c(s);

    auto w = wfc(Radius,0);
    auto c = convol_c2r(sc, w);

    auto Nprj = project_value(c, p0);
    auto Ncic = count_in_sphere(Radius, p, p0);

    for(int i = 0; i < NUMRAN; ++i) ofsprj << Nprj[i] << ", ";
    for(int i = 0; i < NUMRAN; ++i) ofscic << Ncic[i] << ", ";

    const double r0 {0.1};
    const double r1 {100};
    const int NR {30};
    const double dr {(r1 - r0) / NR};

    std::vector<double> rho_p(NR), rho_c(NR);
    std::vector<double> r; for(int i = 0; i < NR; ++i) r[i] = r0 + i * dr;

    
    
    for(int i = 0; i < NR; ++i)
    {
        auto w = wfc(r[i],0);
        auto c = convol_c2r(sc,w);
    }
}

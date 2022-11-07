#include"mracs.h"

#define NUMRAN  100

int main()
{
    read_parameter();
    //auto p1 = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto p2 = read_in_Halo_4vector("/data0/BigMDPL/BigMDPL_halo.bin");

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; for(size_t i = 0; i < NUMRAN; ++i) p0.push_back({u(e), u(e), u(e), 1.});

    force_kernel_type(1);
    //auto cic_dm = count_in_sphere(Radius,p1,p0);
    //auto s1 = sfc(p1);
    auto w = wfc(Radius,0);
    //auto c1 = convol3d(s1,w);
    //auto prj_dm = project_value(c1,p0);
    

    //delete[] s1;
    //delete[] c1;

    auto s2 = sfc(p2);
    auto c2= convol3d(s2,w);
    auto prj_halo = project_value(c2,p0);
    //auto cic_halo = count_in_sphere(Radius,p1,p0);

    //double rho_dm   = p1.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3);
    double rho_halo = p2.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3);
    double eva_dm{0}, eva_halo{0};
    //for(size_t i = 0; i < NUMRAN; ++i) {eva_dm += prj_dm[i] / rho_dm - 1; eva_halo += prj_halo[i] / rho_halo - 1;}
    for(size_t i = 0; i < NUMRAN; ++i) {eva_halo += prj_halo[i] / rho_halo - 1;}
    std::cout << eva_halo << ", " << eva_dm << std::endl;
    std::cout << rho_halo << std::endl;

    //std::cout << "-----------------prj----------------" << std::endl << "===>dm: " << std::endl;
    //for(size_t i = 0; i < NUMRAN; ++i) std::cout << prj_dm[i] / rho_dm - 1 << ", "; std::cout << std::endl << "===>halo: " << std::endl;
    for(size_t i = 0; i < NUMRAN; ++i) std::cout << prj_halo[i]  << ", "; std::cout << std::endl;
    //std::cout << "-----------------cic----------------" << std::endl << "===>dm: " << std::endl;
    //for(size_t i = 0; i < NUMRAN; ++i) std::cout << cic_dm[i] / rho_dm - 1 << ", "; std::cout << std::endl << "===>halo: " << std::endl;
    //for(size_t i = 0; i < NUMRAN; ++i) std::cout << cic_halo[i] / rho_halo - 1 << ", "; std::cout << std::endl;
}
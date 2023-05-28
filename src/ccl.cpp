// Cross-correlation function check
#include"mracs.h"

int main()
{
    read_parameter();
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub005.bin");
    auto p20 = read_in_Halo_4vector("/data0/MDPL2/halo_position.bin");

    std::vector<Particle> p2;  
    const double M_min {2e12};
    for(size_t i = 0; i < p20.size(); ++i) if(p20[i].weight > M_min) p2.push_back({p20[i].x, p20[i].y, p20[i].z, 1.});
    std::vector<Particle>().swap(p20);

    std::cout << "dm: "   << p1.size() << std::endl;
    std::cout << "halo: " << p2.size() << std::endl;
    //#################################################

    auto sc1 = sfc_r2c(sfc(p1),true);
    auto sc2 = sfc_r2c(sfc(p2),true);

    auto ccf  = densityCorrelationDWT(sc1, sc2);
    auto acf1 = densityCorrelationDWT(sc1, sc1);
    auto acf2 = densityCorrelationDWT(sc2, sc2);

    force_base_type(0,1);
    auto ccf_t  = densityCorrelationFFT(sc1, sc2);
    auto acf1_t = densityCorrelationFFT(sc1, sc1);
    auto acf2_t = densityCorrelationFFT(sc2, sc2);


    const double k1 = M_PI * 2 / SimBoxL;
    const int klen = GridLen / 2 + 1;
    std::cout << "k: " << std::endl; for(int i = 0; i < klen; ++i) std::cout << i*k1 << ", ";
    std::cout << std::endl;
    std::cout << "ccf: " << std::endl; for(int i = 0; i < klen; ++i) std::cout << ccf[i]/sqrt(acf1[i]*acf2[i]) << ", ";
    std::cout << std::endl;
    std::cout << "ccft: " << std::endl; for(int i = 0; i < klen; ++i) std::cout << ccf_t[i]/sqrt(acf1_t[i]*acf2_t[i]) << ", ";
    std::cout << std::endl;

}
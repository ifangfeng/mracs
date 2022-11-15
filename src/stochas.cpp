#include"mracs.h"

#define NUMRAN  1000

int main()
{
    read_parameter();
    auto p1 = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto p20 = read_in_Halo_4vector("/data0/BigMDPL/BigMDPL_halo.bin");
    
    std::vector<Particle> p2;  

    const double M_min {2e13};
    for(size_t i = 0; i < p20.size(); ++i) if(p20[i].weight > M_min) p2.push_back({p20[i].x, p20[i].y, p20[i].z, 1.});
    std::vector<Particle>().swap(p20);
    std::cout << p2.size() << std::endl;


    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; for(size_t i = 0; i < NUMRAN; ++i) p0.push_back({u(e), u(e), u(e), 1.});

    double exp_m = p1.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3);
    double exp_h = p2.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3);

    force_kernel_type(1);
    auto w = wfc(Radius,0);

    auto s1 = sfc(p1);
    auto c1 = convol3d(s1,w);
    auto dtm = project_value(c1,p0);
    delete[] s1;
    delete[] c1;

    auto s2 = sfc(p2);
    auto c2= convol3d(s2,w);
    auto dth = project_value(c2,p0);
    delete[] s2;
    delete[] c2;
    delete[] w;

    for(size_t i = 0; i < NUMRAN; ++i){
        dtm[i] = dtm[i] / exp_m - 1;
        dth[i] = dth[i] / exp_h - 1;}
    
    const int num_bin {10};
    const double dtm0 {-1};
    const double dtm1 {4};
    const double ddt = (dtm1 - dtm0) / num_bin;

    for(size_t i = 0; i < NUMRAN; ++i){
        int temp = floor((dtm[i] - dtm0) / ddt);
        if(temp < num_bin && temp >= 0)

    }

    
    double ave_dm{0}, ave_halo{0};
    for(size_t i = 0; i < NUMRAN; ++i) {ave_dm += dtm[i]; ave_halo += dth[i];}
    std::cout << ave_halo / NUMRAN << ", " << ave_dm / NUMRAN << std::endl;
    std::cout << exp_h << ", " << exp_m << std::endl;

    std::cout << "-----------------prj----------------" << std::endl << "===>dm: " << std::endl;
    for(size_t i = 0; i < NUMRAN; ++i) std::cout << dtm[i] / exp_m - 1 << ", "; std::cout << std::endl << "===>halo: " << std::endl;
    for(size_t i = 0; i < NUMRAN; ++i) std::cout << dth[i] / exp_h - 1 << ", "; std::cout << std::endl;





    /*for(size_t i = 0; i < NUMRAN; ++i) 
    {
        if(prj_halo[i] < -exp_halo)
            std::cout << "halo"<< p0[i].x << ",  "<< p0[i].y << ", "<< p0[i].z << ",  |"<< prj_halo[i]<< ", "<<std::endl;
    };
    for(size_t i = 0; i < NUMRAN; ++i) 
    {
        if(prj_dm[i] < -exp_dm)
            std::cout << "dm" << p0[i].x << ",  "<< p0[i].y << ", "<< p0[i].z << ",  |"<< prj_dm[i]<< ", "<<std::endl;
    };
    */

    //auto cic_dm = count_in_sphere(Radius,p1,p0);
    //auto cic_halo = count_in_sphere(Radius,p1,p0);
    //std::cout << "-----------------cic----------------" << std::endl << "===>dm: " << std::endl;
    //for(size_t i = 0; i < NUMRAN; ++i) std::cout << cic_dm[i] / exp_dm - 1 << ", "; std::cout << std::endl << "===>halo: " << std::endl;
    //for(size_t i = 0; i < NUMRAN; ++i) std::cout << cic_halo[i] / exp_halo - 1 << ", "; std::cout << std::endl;
}

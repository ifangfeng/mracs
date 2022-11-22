#include"mracs.h"

#define NUMRAN  1000

int main()
{
    std::string ofname_dm_cic {"output/dm1.txt"};
    std::string ofname_halo_cic {"output/halo1.txt"};
    std::string ofname_dm_prj {"output/dm0.txt"};
    std::string ofname_halo_prj {"output/halo0.txt"};

    std::ofstream ofsm1 (ofname_dm_cic);
    std::ofstream ofsh1 (ofname_halo_cic);
    std::ofstream ofsm0 (ofname_dm_prj);
    std::ofstream ofsh0 (ofname_halo_prj);

    if(!ofsm1 || !ofsh1 || !ofsm0 || !ofsh0)
    {
        std::cout << "Abort! opening write out file with error..." << std::endl;
        std::terminate();
    }

    read_parameter();
    auto p10 = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto p20 = read_in_Halo_4vector("/data0/BigMDPL/BigMDPL_halo.bin");
    
    std::vector<Particle> p1, p2;  

    for(size_t i = 0; i < p10.size(); i += 1) p1.push_back({p10[i].x, p10[i].y, p10[i].z, 1.});
    std::vector<Particle>().swap(p10);
    std::cout << "dm: " << p1.size() << std::endl;

    const double M_min {2e12};
    for(size_t i = 0; i < p20.size(); ++i) if(p20[i].weight > M_min) p2.push_back({p20[i].x, p20[i].y, p20[i].z, 1.});
    std::vector<Particle>().swap(p20);
    std::cout << "halo: " << p2.size() << std::endl;
    //#################################################


    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; for(size_t i = 0; i < NUMRAN; ++i) p0.push_back({u(e), u(e), u(e), 1.});

    double exp_m = p1.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3); std::cout << "exp_m: " << exp_m << std::endl;
    double exp_h = p2.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3); std::cout << "exp_h: " << exp_h << std::endl;

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
    delete[] w;
    delete[] s2;
    delete[] c2;
    
    auto dtmcic = count_in_sphere(Radius, p1, p0);
    auto dthcic = count_in_sphere(Radius, p2, p0);


    for(size_t i = 0; i < NUMRAN; ++i){
        dtm[i]    = dtm[i] / exp_m - 1;
        dth[i]    = dth[i] / exp_h - 1;
        dtmcic[i] = dtmcic[i] / exp_m - 1;
        dthcic[i] = dthcic[i] / exp_h - 1;
    }

    
    
    const int num_bin {15};
    const double dtm0 {-1}, dtm1 {2};
    const double ddt {(dtm1 - dtm0) / num_bin};

    std::vector<unsigned> count(num_bin);
    std::vector<double> ave(num_bin), var(num_bin), cbin(num_bin);
    for(int i = 0; i < num_bin; ++i) {ave[i] = 0; var[i] = 0; count[i] = 0;}
    for(int i = 0; i < num_bin; ++i) {cbin[i] = dtm0 + (i + 0.5) * ddt;}

    std::vector<unsigned> countcic(num_bin);
    std::vector<double> avecic(num_bin), varcic(num_bin);
    for(int i = 0; i < num_bin; ++i) {avecic[i] = 0; varcic[i] = 0; countcic[i] = 0;}

    for(size_t i = 0; i < NUMRAN; ++i){
        int index = floor((dtm[i] - dtm0) / ddt);
        if(index < num_bin && index >= 0){
            ave[index] += dth[i];
            var[index] += pow(dth[i],2);
            ++count[index];
        }
    }

    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]) {
            ave[i] /= count[i];
            var[i] /= count[i];
            var[i] -= pow(ave[i],2);
        }
    }

    for(size_t i = 0; i < NUMRAN; ++i){
        int index = floor((dtmcic[i] - dtm0) / ddt);
        if(index < num_bin && index >= 0){
            avecic[index] += dthcic[i];
            varcic[index] += pow(dthcic[i],2);
            ++countcic[index];
        }
    }

    for(size_t i = 0; i < num_bin; ++i){
        if(countcic[i]) {
            avecic[i] /= countcic[i];
            varcic[i] /= countcic[i];
            varcic[i] -= pow(avecic[i],2);
        }
    }

    std::cout << "centre of bin: " << std::endl; for(auto i : cbin) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "=============prj=============";
    std::cout << "count in bin: " << std::endl; for(auto i : count) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h: " << std::endl; for(auto i : ave) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h deviation: " << std::endl; for(auto i : var) std::cout << sqrt(i) << ", "; std::cout << std::endl;
    std::cout << "=============cic=============";
    std::cout << "count in bin: " << std::endl; for(auto i : countcic) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h: " << std::endl; for(auto i : avecic) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h deviation: " << std::endl; for(auto i : varcic) std::cout << sqrt(i) << ", "; std::cout << std::endl;
    
    for(size_t i = 0; i < NUMRAN; ++i) ofsm0 << dtm[i] << ", ";
    for(size_t i = 0; i < NUMRAN; ++i) ofsh0 << dth[i] << ", ";
    for(size_t i = 0; i < NUMRAN; ++i) ofsm1 << dtmcic[i] << ", ";
    for(size_t i = 0; i < NUMRAN; ++i) ofsh1 << dthcic[i] << ", ";

    std::cout << std::endl;
    
}

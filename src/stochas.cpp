#include"mracs.h"

#define NUMRAN  10000

int main()
{
    std::string ofname_dm {"dm.txt"};
    std::string ofname_halo {"halo.txt"};
    std::string ofname_header {"header.txt"};

    std::ofstream ofsm (ofname_dm);
    std::ofstream ofsh (ofname_halo);
    std::ofstream ofsd (ofname_header);

    if(!ofsm || !ofsh || !ofsd)
    {
        std::cout << "Abort! opening write out file with error..." << std::endl;
        std::terminate();
    }

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

    //for(size_t i = 0; i < 1e6; ++i) {p1.push_back({u(e), u(e), u(e), 1.}); p2.push_back({u(e), u(e), u(e), 1.});}

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

    std::cout << "hello1" << std::endl;

    for(size_t i = 0; i < NUMRAN; ++i){
        dtm[i] = dtm[i] / exp_m - 1;
        dth[i] = dth[i] / exp_h - 1;}
    
    const int num_bin {15};
    const double dtm0 {-0.5}, dtm1 {0.5};
    const double ddt {(dtm1 - dtm0) / num_bin};

    std::vector<unsigned> count(num_bin);
    std::vector<double> ave(num_bin), var(num_bin), cbin(num_bin);
    for(int i = 0; i < num_bin; ++i) {ave[i] = 0; var[i] = 0; count[i] = 0;}
    for(int i = 0; i < num_bin; ++i) {cbin[i] = dtm0 + (i + 0.5) * ddt;}


    std::cout << "hello2" << std::endl;

    for(size_t i = 0; i < NUMRAN; ++i){
        int index = floor((dtm[i] - dtm0) / ddt);
        if(index < num_bin && index >= 0){
            ave[index] += dth[i];
            var[index] += pow(dth[i],2);
            ++count[index];
        }
    }
    std::cout << "hello3" << std::endl;
    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]) {
            ave[i] /= count[i];
            var[i] /= count[i];
            var[i] -= pow(ave[i],2);
        }
    }

    std::cout << "count in bin: " << std::endl; for(auto i : count) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "centre of bin: " << std::endl; for(auto i : cbin) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h: " << std::endl; for(auto i : ave) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h deviation: " << std::endl; for(auto i : var) std::cout << sqrt(i) << ", "; std::cout << std::endl;
    std::cout << "writing text file...";
    
    for(size_t i = 0; i < NUMRAN; ++i) ofsm << dtm[i] << ", "; ofsm.close();
    for(size_t i = 0; i < NUMRAN; ++i) ofsh << dth[i] << ", "; ofsh.close();
    for(auto i : count) ofsd << i << ", "; ofsd << std::endl;
    for(auto i : cbin) ofsd << i << ", "; ofsd << std::endl;
    for(auto i : ave) ofsd << i << ", "; ofsd << std::endl;
    for(auto i : var) ofsd << sqrt(i) << ", "; ofsd.close();
    std::cout << std::endl;
    
}

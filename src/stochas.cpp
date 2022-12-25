#include"mracs.h"

#define NUMRAN  5000
void do_something
(std::string cic_or_prj, std::ofstream& ofsm, std::ofstream& ofsh, double* dtm, double* dth, double expect_m, double expect_h, double R);

int main()
{
    read_parameter();
    force_kernel_type(1);
    const double M_min {2e12};
    auto p1  = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2_raw = read_in_Halo_4vector("/data0/MDPL2/halo_position.bin");

    std::vector<Particle> p2;for(auto i : p2_raw) if(i.weight > M_min) p2.push_back({i.x, i.y, i.z, 1.});
    std::vector<Particle>().swap(p2_raw);
    std::cout << "dm: "   << p1.size() << std::endl;
    std::cout << "halo: " << p2.size() << std::endl;
    //#################################################

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(51, SimBoxL-51);
    std::vector<Particle> p0; for(size_t i = 0; i < NUMRAN; ++i) p0.push_back({u(e), u(e), u(e), 1.});

    //#################################################
    std::vector<int> R_vec {15};//50,30,20,15,10,5}; 
    std::vector<int> J_vec {9};//{8,7,6,5,4};
    for(auto R : R_vec){
        double exp_m = p1.size() * 4./3 * M_PI * pow(R / SimBoxL,3); 
        double exp_h = p2.size() * 4./3 * M_PI * pow(R / SimBoxL,3); 
        std::cout << "exp_m: " << exp_m << std::endl;
        std::cout << "exp_h: " << exp_h << std::endl;

        for(auto J : J_vec){
            std::string ofname_prj_dm = "output/prj_R" + std::to_string(R) + "J" + std::to_string(J) + "_dm.txt";
            std::string ofname_prj_halo  = "output/prj_R" + std::to_string(R) + "J" + std::to_string(J) + "_halo.txt";
            std::ofstream ofs_prj_dm {ofname_prj_dm};
            std::ofstream ofs_prj_halo {ofname_prj_halo};

            force_resoluton_J(J);
            auto w = wfc(R, 0);
            auto s1 = sfc(p1); 
            auto s2 = sfc(p2); 
            auto c1 = convol3d(s1,w); 
            auto c2 = convol3d(s2,w);
            auto dtmprj = project_value(c1,p0); delete[] s1; delete[] c1;
            auto dthprj = project_value(c2,p0); delete[] s2; delete[] c2; delete[] w;

            std::cout << "================================prj_R" << R << "J" << J << "================: " << std::endl;
            do_something("prj", ofs_prj_dm, ofs_prj_halo, dtmprj, dthprj, exp_m, exp_h, R);
        }
    }
    for(auto R : R_vec){
        std::string ofname_cic_dm   = "output/cic_R" + std::to_string(R) + "_dm.txt"; 
        std::string ofname_cic_halo = "output/cic_R" + std::to_string(R) + "_halo.txt"; 
        std::ofstream ofs_cic_dm   {ofname_cic_dm};
        std::ofstream ofs_cic_halo {ofname_cic_halo};
        
        double exp_m = p1.size() * 4./3 * M_PI * pow(R / SimBoxL,3); 
        double exp_h = p2.size() * 4./3 * M_PI * pow(R / SimBoxL,3); 
        std::cout << "exp_m: " << exp_m << std::endl;
        std::cout << "exp_h: " << exp_h << std::endl;
    
        auto dtmcic = count_in_sphere(R, p1, p0);
        auto dthcic = count_in_sphere(R, p2, p0);
        
        std::cout << "++++++++++++++++++++++++++++++++++++cic_R" << R << "++++++++++++++++: " << std::endl;
        do_something("cic", ofs_cic_dm, ofs_cic_halo, dtmcic,dthcic, exp_m, exp_h, R);
    }
}


void do_something
(std::string cic_or_prj, std::ofstream& ofsm, std::ofstream& ofsh, double* dtm, double* dth, double expect_m, double expect_h, double R){
    for(size_t i = 0; i < NUMRAN; ++i){
        dtm[i] = dtm[i]/expect_m - 1;
        dth[i] = dth[i]/expect_h - 1;
    }

    const int num_bin {15};
    double dtm0 {-1}, dtm1 {2};
    if(R < 20){
        dtm0 = -1;
        dtm1 = 2;
    }
    else if(R < 30){
        dtm0 = -1;
        dtm1 = 1.5;
    }
    else if(R < 50){
        dtm0 = -0.5;
        dtm1 = 1;
    }
    else if(R <= 80){
        dtm0 = -0.4;
        dtm1 = 0.5;
    }
    

    const double ddt {(dtm1 - dtm0) / num_bin};

    std::vector<unsigned> count(num_bin);
    std::vector<double> ave(num_bin), var(num_bin), cbin(num_bin);
    for(int i = 0; i < num_bin; ++i) {ave[i] = 0; var[i] = 0; count[i] = 0;}
    for(int i = 0; i < num_bin; ++i) {cbin[i] = dtm0 + (i + 0.5) * ddt;}

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
    double x_tmp{0}, y_tmp{0}, bias1,bias2;
    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]){
            x_tmp += cbin[i];
            y_tmp += ave[i];
        }
    }
    bias1 = y_tmp/x_tmp;
    x_tmp = 0;
    y_tmp = 0;
    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]){
            x_tmp += cbin[i] * count[i];
            y_tmp += ave[i] * count[i];
        }
    }
    bias2 = y_tmp/x_tmp;
    std::cout << "bias1= " << bias1 << std::endl;
    std::cout << "bias2= " << bias2 << std::endl;
    std::cout << "centre of bin: " << std::endl; for(auto i : cbin) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl; for(auto i : count) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h: " << std::endl; for(auto i : ave) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h deviation: " << std::endl; for(auto i : var) std::cout << sqrt(i) << ", "; std::cout << std::endl;

    for(size_t i = 0; i < NUMRAN; ++i) ofsm << dtm[i] << ", ";
    for(size_t i = 0; i < NUMRAN; ++i) ofsh << dth[i] << ", ";

    delete[] dtm;
    delete[] dth;
    std::cout << std::endl;
}


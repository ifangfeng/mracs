// Probability Distribution Function of [density filed]
#include"mracs.h"

//#define halo
void pdf(std::vector<Particle> p, std::vector<Particle> p1);
int64_t npartall{0};

int main()
{
    read_parameter();
    #ifdef halo
    auto p10 = read_in_Halo_4vector("/data0/MDPL2/halo_position.bin");
    std::vector<Particle> p1;  
    const double M_min {2e12};
    for(size_t i = 0; i < p10.size(); ++i) if(p10[i].weight > M_min) p1.push_back({p10[i].x, p10[i].y, p10[i].z, 1.});
    std::vector<Particle>().swap(p10);
    std::cout << "halo: " << p1.size() << std::endl;
    #else 
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub005.bin");
    #endif

    npartall = p1.size();

    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int rand = 1000;
    for(int i = 0; i < rand; ++i)
        for(int j = 0; j < rand; ++j)
            for(int k = 0; k < rand; ++k)
                p.push_back({(i + u(e)) / rand * SimBoxL, (j + u(e)) / rand * SimBoxL, (k + u(e)) / rand * SimBoxL});

    std::vector<int> resolvec {7,8,9,10,11};
    for(auto i : resolvec){
        force_resoluton_J(i);
        pdf(p,p1);
    }
}


void pdf(std::vector<Particle> p, std::vector<Particle> p1)
{
    auto s = sfc(p1);
    auto w = wfc(Radius,0);
    auto c = convol3d(s,w);
    delete[] w;
    delete[] s;


    const int nbin {500};
    const double rhomax {10};
    const double rhomin {0};
    const double d_rho = (rhomax - rhomin) / nbin;
    const double cicexpect = npartall * 4./3 * M_PI * pow(Radius/SimBoxL, 3);
    std::cout << "cic expectation: " << cicexpect << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = (i + 0.5) * d_rho + rhomin;
    double count[nbin]{0};
    double value[nbin]{0};
    double sigma[nbin]{0};
    auto n_prj = project_value(c,p);
    for(size_t i = 0; i < p.size(); ++i) {
        int index = (n_prj[i] / cicexpect - rhomin) / d_rho;
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / p.size() / d_rho;
    }
    /*
    const int NJK{100}; // number of Jack knife subsample
    const uint64_t sublen {p.size() / NJK};
    double c_tmp[nbin]{0};
    double c_ave[nbin]{0};
    for(int i = 0; i < NJK; ++i){
        for(int n = 0; n < sublen; ++n){
            int index = n_prj[i + n * NJK] / cicexpect / d_rho;
            if(index < nbin && index >= 0)
            ++c_tmp[index];
        }
        for(int j = 0; j < nbin; ++j){
            c_tmp[j] /= sublen * d_rho;
            sigma[j] += pow(c_tmp[j], 2);
            c_ave[j] += c_tmp[j];
            c_tmp[j] = 0;
        }
    }
    for(int i = 0; i < nbin; ++i){
        sigma[i] = sqrt((sigma[i] - c_ave[i]) / (NJK - 1));
    }
    */
    std::cout << "========out put========J" << Resolution << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl;
    //std::cout << "standard deviation: " << '\n';for(int i = 0; i < nbin; ++i) std::cout << sigma[i] << ", "; std::cout << std::endl;
}
// Probability Distribution Function of [density filed]
#include"mracs.h"

//#define halo
void pdf(double* c, size_t N);
void cic_pdf(double* c, size_t N);
void cic_pdf_all(double* c, size_t N);
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
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    #endif

    npartall = p1.size();

    std::vector<Particle> p;
    std::default_random_engine e;
    /*
    std::uniform_real_distribution<double> u(0,1);
    const int rand = 1000;
    for(int i = 0; i < rand; ++i)
        for(int j = 0; j < rand; ++j)
            for(int k = 0; k < rand; ++k)
                p.push_back({(i + u(e)) / rand * SimBoxL, (j + u(e)) / rand * SimBoxL, (k + u(e)) / rand * SimBoxL});
    */
    std::uniform_real_distribution<double> u(50.1, SimBoxL-50.1);
    for(int i = 0; i < 1000; ++i) p.push_back({u(e),u(e),u(e),1.});

    std::vector<int> resolvec {5,6,7,8,9};
    for(auto i : resolvec){
        force_resoluton_J(i);
        auto s = sfc(p1);
        auto w = wfc(Radius,0);
        auto c = convol3d(s,w);
        delete[] w;
        delete[] s;
        auto n_prj = project_value(c,p);
        pdf(n_prj, p.size());
    }

    auto n_cic = count_in_sphere(Radius,p1,p);
    pdf(n_cic, p.size());


    
}


void pdf(double* c, size_t N)
{
    const int nbin {168};
    const double rhomax {5};
    const double rhomin {0};
    const double d_rho = (rhomax - rhomin) / nbin;
    const double cicexpect = npartall * 4./3 * M_PI * pow(Radius/SimBoxL, 3);
    std::cout << "cic expectation: " << cicexpect << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = (i + 0.5) * d_rho + rhomin;
    double count[nbin]{0};
    double value[nbin]{0};
    for(size_t i = 0; i < N; ++i) {
        int index = (c[i] / cicexpect - rhomin) / d_rho;
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / N / d_rho;
    }
    
    const int NJK{100}; // number of Jack knife subsample
    const uint64_t sublen {N / NJK};
    double count_tmp[nbin]{0};
    double value_tmp[nbin]{0};
    double xsquare[nbin]{0};
    double xbar[nbin]{0};
    double sigma[nbin]{0};
    double sigma_percent[nbin]{0};
    std::cout << sublen << std::endl;
    for(int n = 0; n < NJK; ++n)
    {
        #pragma omp paralle for reduction(+:count_tmp)
        for(int i = 0; i < N; ++i)
        {
            if(i%NJK != n)
            {
                int index = (c[i] / cicexpect  - rhomin) / d_rho;
                if(index < nbin && index >= 0)
                ++count_tmp[index];
            }
        }
        for(int j = 0; j < nbin; ++j){
            value_tmp[j] = count_tmp[j] / sublen / d_rho;
            xsquare[j] += pow(value_tmp[j], 2);
            xbar[j] += value_tmp[j];
            count_tmp[j] = 0;
        }
    }
    for(int i = 0; i < nbin; ++i){
        sigma[i] = sqrt((xsquare[i] - pow(xbar[i],2)/NJK) / (NJK - 1));
        if(value[i])
        sigma_percent[i] = 100 * sigma[i] / value[i];
    }
    
    std::cout << "========out put========J" << Resolution << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl;
    std::cout << "sigma:: " << std::endl;       for(int i = 0; i < nbin; ++i) std::cout << sigma[i] << ", "; std::cout << std::endl;
    std::cout << "sigma percent: " << std::endl;for(int i = 0; i < nbin; ++i) std::cout << sigma_percent[i] << ", "; std::cout << std::endl;
}

void cic_pdf(double* c, size_t N)
{
    const double rhomax {5};
    const double rhomin {0};
    const double cicexpect = npartall * 4./3 * M_PI * pow(Radius/SimBoxL, 3);
    const double d_rho = 1./cicexpect;
    const int nbin = (rhomax - rhomin) / d_rho + 1;

    std::cout << "cic expectation: " << cicexpect << std::endl;
    std::cout << "number of bin: " << nbin << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = i * d_rho + rhomin;
    double count[nbin]{0};
    double value[nbin]{0};
    for(size_t i = 0; i < N; ++i) {
        int index = static_cast<int>(c[i]);
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / N / d_rho;
    }
    
    const int NJK{100}; // number of Jack knife subsample
    const uint64_t sublen {N / NJK};
    double count_tmp[nbin]{0};
    double value_tmp[nbin]{0};
    double xsquare[nbin]{0};
    double xbar[nbin]{0};
    double sigma[nbin]{0};
    double sigma_percent[nbin]{0};
    std::cout << sublen << std::endl;
    for(int n = 0; n < NJK; ++n)
    {
        #pragma omp paralle for reduction(+:count_tmp)
        for(int i = 0; i < N; ++i)
        {
            if(i%NJK != n)
            {
                int index = static_cast<int>(c[i]);
                if(index < nbin && index >= 0)
                ++count_tmp[index];
            }
        }
        for(int j = 0; j < nbin; ++j){
            value_tmp[j] = count_tmp[j] / sublen / d_rho;
            xsquare[j] += pow(value_tmp[j], 2);
            xbar[j] += value_tmp[j];
            count_tmp[j] = 0;
        }
    }
    for(int i = 0; i < nbin; ++i){
        sigma[i] = sqrt((xsquare[i] - pow(xbar[i],2)/NJK) / (NJK - 1));
        if(value[i])
        sigma_percent[i] = 100 * sigma[i] / value[i];
    }
    
    std::cout << "========out put========cic" << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl;
    std::cout << "sigma:: " << std::endl;       for(int i = 0; i < nbin; ++i) std::cout << sigma[i] << ", "; std::cout << std::endl;
    std::cout << "sigma percent: " << std::endl;for(int i = 0; i < nbin; ++i) std::cout << sigma_percent[i] << ", "; std::cout << std::endl;
}






void cic_pdf_all(double* c, size_t N){
    int count [200]{0};
    double cbin [200];
    const double cicexpect = npartall * 4./3 * M_PI * pow(Radius/SimBoxL, 3);
    for(int i = 0; i < 200; ++i)
    {
        cbin[i] = i/cicexpect;
    }
    for(int n = 0; n < N; ++n)
    {
        ++count[static_cast<int>(c[n])];
    }
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < 200; ++i) std::cout << cbin[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < 200; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    

}
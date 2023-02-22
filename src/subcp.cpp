// Probability Distribution Function of [density filed]
#include"mracs.h"
#include"kdtree.hpp"

void write_1MP_to_file(double* ncount, std::ofstream& ofs);
void pdf(double* c, size_t N, size_t npartall);
void cic_pdf(double* c, size_t N, size_t npartall);
void cic_pdf_all(double* c, size_t N, size_t npartall);

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
    auto p1 = read_in_DM_3vector(DataDirec);
    
    std::vector<Particle> p0;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int nps = 100;
    const double safeband = 50;

    for(int i = 0; i < nps; ++i)
        for(int j = 0; j < nps; ++j)
            for(int k = 0; k < nps; ++k)
                p0.push_back({safeband + (i + u(e)) * (SimBoxL - 2*safeband) / nps,
                              safeband + (j + u(e)) * (SimBoxL - 2*safeband) / nps,
                              safeband + (k + u(e)) * (SimBoxL - 2*safeband) / nps});

    auto s = sfc(p1);
    auto w = wfc(Radius,0);
    auto c = convol3d(s,w);
    delete[] w;
    delete[] s;
    auto n_prj = project_value(c,p0);
    write_1MP_to_file(n_prj, ofsprj);
    pdf(n_prj, p0.size(), p1.size());

    // ==================================================
    // construction of kd-tree
    Kdtree::KdNodeVector nodes;
    double diff;
    clock_t begin, end;
    for (size_t i = 0; i < p1.size(); ++i) {
      std::vector<double> point(3);
      point[0] = p1[i].x;
      point[1] = p1[i].y;
      point[2] = p1[i].z;
      nodes.push_back(Kdtree::KdNode(point));
    }
    begin = clock();
    Kdtree::KdTree tree(&nodes);
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Creation time for " << p1.size() << " galaxy points:  " << diff << "s" << std::endl;

    // 1.2) range query
    Kdtree::KdNodeVector result;
    auto n_cic = new double[p0.size()];
    begin = clock();
    for (size_t i = 0; i < p0.size(); ++i) {
      std::vector<double> test_point(3);
      test_point[0] = p0[i].x;
      test_point[1] = p0[i].y;
      test_point[2] = p0[i].z;
      tree.range_nearest_neighbors(test_point, Radius, &result);
      n_cic[i] = result.size();
    }
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "counting time(tree) for " << p0.size() << " random points:  " << diff << "s" << std::endl;
    write_1MP_to_file(n_prj, ofsprj);
    pdf(n_cic, p0.size(), p1.size());


    
}

void write_1MP_to_file(double* ncount, std::ofstream& ofs)
{
    double tmp{0};
    for(int i = 0; i < 1000000; ++i){
        ofs << ncount[i];
    }
}

void pdf(double* c, size_t N, size_t npartall)
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

void cic_pdf(double* c, size_t N, size_t npartall)
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






void cic_pdf_all(double* c, size_t N, size_t npartall){
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
#include"mracs.h"


std::vector<Particle> snap_push_back(std::string ifname);
void sub_all(std::vector<Particle> p0);
void sub_particle(std::vector<Particle> p, std::vector<Particle> ps);
void array_merge(double* a, double* b, size_t N);
double* init_array(size_t N);
void load_and_clean(std::string ifname, double* s0);

int64_t npartall{0};
std::vector<Particle> p1,p2,p3,p4;


int main()
{
    read_parameter();
    std::string direct1 {"/mnt/disk/wd/MDPL2/"};
    std::string direct2 {"/mnt/disk/data/MDPL2/snap_130/"};
    std::string fname;
    auto s0 = init_array(GridNum);
    for(int n = 0; n < 20; ++n)
    {
        if(n < 11){
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 100; i += 1){
                fname = direct1 + "snap_130." + std::to_string(100*n+i);
                load_and_clean(fname, s0);}
        }
        else if (n < 19){
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 100; i += 1){
                fname = direct2 + "snap_130." + std::to_string(100*n+i);
                load_and_clean(fname, s0);}
        }
        else{
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 20; i += 1){
                fname = direct2 + "snap_130." + std::to_string(100*n+i);
                load_and_clean(fname, s0);}
        }
    }
    std::cout << "-----> Total number of particles: " << npartall << " = " << pow(npartall, 1./3) << "^3" << std::endl;
    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    for(int i = 0; i < 100; ++i)
        for(int j = 0; j < 100; ++j)
            for(int k = 0; k < 100; ++k)
                p.push_back({i + u(e),j + u(e), k + u(e)});

    auto w  = wfc(Radius,0);
    auto s1 = sfc(p1);        std::vector<Particle>().swap(p1);
    auto s2 = sfc(p2);        std::vector<Particle>().swap(p2);
    auto s3 = sfc(p3);        std::vector<Particle>().swap(p3);
    auto s4 = sfc(p4);        std::vector<Particle>().swap(p4);
    auto c0 = convol3d(s0,w); delete[] s0;
    auto c1 = convol3d(s1,w); delete[] s1;
    auto c2 = convol3d(s2,w); delete[] s2;
    auto c3 = convol3d(s3,w); delete[] s3;
    auto c4 = convol3d(s4,w); delete[] s4;
    delete[] w;

    const int nbin {20};
    const double rhomax {4};
    const double d_rho = rhomax / nbin;
    const double cicexpect0 = npartall * 4./3 * M_PI * pow(Radius/SimBoxL, 3);
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = (i + 0.5) * d_rho;
    double count[nbin]{0};
    double value[nbin]{0};
    double sigma[nbin]{0};
    auto n_prj0 = project_value(c0,p);
    for(size_t i = 0; i < p.size(); ++i) 
    {
        int index = n_prj0[i] / cicexpect0 / d_rho;
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / p.size() / d_rho;
    }
    const int NJK{100}; // number of Jack knife subsample
    double c_tmp[nbin]{0};
    double c_ave[nbin]{0};
    for(int i = 0; i < NJK; ++i)
    {
        for(int j = 0; j < p.size() / NJK; ++j)
        {
            int index = n_prj0[i + j * NJK] / cicexpect0 / d_rho;
            ++c_tmp[index];
        }
        for(int i = 0; i < nbin; ++i)
        {
            c_tmp[i] /= npartall * d_rho;
            sigma[i] += pow(c_tmp[i],2);
            c_ave[i] += c_tmp[i];
        }
    }
    for(int i = 0; i < nbin; ++i)
    {
        sigma[i] = sqrt((sigma[i] - c_ave[i]) / (NJK - 1));
    }


    auto n_prj1 = project_value(c1,p);
    auto n_prj2 = project_value(c2,p);
    auto n_prj3 = project_value(c3,p);
    auto n_prj4 = project_value(c4,p);

}


void load_and_clean(std::string ifname, double* s0)
{
    auto p_tmp = snap_push_back(ifname);
    sub_all(p_tmp);
    auto s_tmp = sfc(p_tmp);
    array_merge(s0,s_tmp,GridNum);
    std::vector<Particle>().swap(p_tmp);
    delete[] s_tmp;
}

double* init_array(size_t N)
{
    auto a = new double[N];
    #pragma omp parallel for
    for(size_t i = 0; i < N; ++i)
        a[i] = 0;
    return a;
}

void array_merge(double* a, double* b, size_t N)
{
    #pragma omp parallel for
    for(size_t i = 0; i < N; ++i)
        a[i] += b[i];
}

void sub_all(std::vector<Particle> p0)
{
    sub_particle(p0,p1);
    sub_particle(p1,p2);
    sub_particle(p2,p3);
    sub_particle(p3,p4);
}

void sub_particle(std::vector<Particle> p, std::vector<Particle> ps)
{
    for(size_t i = 0; i < p.size(); i += 10)
    {
        ps.push_back({p[i].x, p[i].y, p[i].z, 1.});
    }
}

std::vector<Particle> snap_push_back(std::string ifname)
{
    std::vector<Particle> p;
    std::ifstream ifs(ifname, std::ios_base::binary);
    ifs.seekg(4 + 256 + 4);
    int nbyte;
    int npart;
    ifs.read(as_bytes(nbyte),sizeof(nbyte));
    npart = nbyte / 4 / 3;
    npartall += npart;
    float a[3];
    for(int i = 0; i < npart; ++i){
        ifs.read(as_bytes(a),sizeof(a));
        p.push_back({a[0],a[1],a[2],1.});
    }
    return p;
}

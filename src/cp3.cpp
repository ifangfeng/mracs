#include"mracs.h"


std::vector<Particle> snap_push_back(std::string ifname);
void sub_all(std::vector<Particle> p0);
void sub_particle(std::vector<Particle> p, std::vector<Particle> ps);
void array_merge(double* a, double* b, size_t N);
double* init_array(size_t N);
void load_and_clean(std::string ifname, double* s0);
double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport);

const int phiSupport = phi[phi.size() - 1] - phi[phi.size() - 2];
std::vector<int> step(phiSupport);
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
    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    for(int i = 0; i < 100; ++i)
        for(int j = 0; j < 100; ++j)
            for(int k = 0; k < 100; ++k)
                p.push_back({i + u(e),j + u(e), k + u(e)});
    auto s1 = sfc(p1);
    auto s2 = sfc(p2);
    auto s3 = sfc(p3);
    auto s4 = sfc(p4);
    auto w = wfc(Radius,0);
    auto c0 = convol3d(s0,w);
    auto c1 = convol3d(s1,w);
    auto c2 = convol3d(s2,w);
    auto c3 = convol3d(s3,w);
    auto c4 = convol3d(s4,w);
    delete[] s0;
    delete[] s1;
    delete[] s2;
    delete[] s3;
    delete[] s4;
    delete[] w;

    for(int i = 0; i < phiSupport; ++i) step[i] = i * SampRate;

    const int bin {20};
    double dens[bin];
    double sigma[bin];
    double count_in_bin[bin];
    for(auto i : p){

    }



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
    float a[3];
    for(int i = 0; i < npart; ++i){
        ifs.read(as_bytes(a),sizeof(a));
        p.push_back({a[0],a[1],a[2],1.});
    }
    return p;
}

double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport)
{
    double sum{0};
    int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
    int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
    int zzc = floor(zz), zzf = SampRate * (zz - zzc);
    
    for(int i = 0; i < phiSupport; ++i)
        for(int j = 0; j < phiSupport; ++j)
            for(int k = 0; k < phiSupport; ++k)
            {
                sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                        * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
            }
    return sum;
}
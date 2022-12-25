// Probability Distribution Function of [density filed]
#include"mracs.h"


void pdf(std::vector<Particle>& p, double* s, double R);
int64_t npartall{0};

void snap_push_back(std::string ifname, std::vector<Particle>& p);
void array_merge(double* a, double* b, size_t N);
double* init_array(size_t N);
void merge_and_clean(std::vector<Particle>& p, double* s);

std::vector<Particle> p1,p2,p3;

int main()
{
    read_parameter();
    std::string direct1 {"/mnt/disk/wd/MDPL2/"};
    std::string direct2 {"/mnt/disk/data/MDPL2/snap_130/"};
    std::string fname;

    auto s0 = init_array(GridNum);

    for(int n = 19; n >= 0; --n)
    {
        if(n < 11){
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 100; i += 1){
                fname = direct1 + "snap_130." + std::to_string(100*n+i);
                snap_push_back(fname, p1);
            }
            merge_and_clean(p1, s0);
        }
        else if (n < 19){
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 100; i += 1){
                fname = direct2 + "snap_130." + std::to_string(100*n+i);
                snap_push_back(fname, p2); 
            }
            merge_and_clean(p2, s0);
        }
        else{
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>n = " << n << std::endl;
            for(int i = 0; i < 20; i += 1){
                fname = direct2 + "snap_130." + std::to_string(100*n+i);
                snap_push_back(fname, p3); 
            }
            merge_and_clean(p3, s0);
        }
    }
    std::cout << "-----> Total number of particles: " << npartall << " = " << pow(npartall, 1./3) << "^3" << std::endl;

    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int rand = 1000;
    for(int i = 0; i < rand; ++i)
        for(int j = 0; j < rand; ++j)
            for(int k = 0; k < rand; ++k)
                p.push_back({(i + u(e)) / rand * SimBoxL, (j + u(e)) / rand * SimBoxL, (k + u(e)) / rand * SimBoxL});

    std::vector<double> Rvec {1,3,5};
    for(auto x : Rvec)
        pdf(p,s0,x);
}


void pdf(std::vector<Particle>& p, double* s, double R)
{
    auto w = wfc(R,0);
    auto c = convol3d(s,w);
    delete[] w;

    const int nbin {500};
    const double rhomax {5};
    const double rhomin {0};
    const double d_rho = (rhomax - rhomin) / nbin;
    const double cicexpect = npartall * 4./3 * M_PI * pow(R/SimBoxL, 3);
    std::cout << "cic expectation: " << cicexpect << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = (i + 0.5) * d_rho + rhomin;
    double count[nbin]{0};
    double value[nbin]{0};
    double sigma[nbin]{0};
    auto n_prj = project_value(c,p);
    delete[] c;
    for(size_t i = 0; i < p.size(); ++i) {
        int index = (n_prj[i] / cicexpect - rhomin) / d_rho;
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / p.size() / d_rho;
    }

    std::cout << "========out put========R" << R << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl;
    //std::cout << "standard deviation: " << '\n';for(int i = 0; i < nbin; ++i) std::cout << sigma[i] << ", "; std::cout << std::endl;
}

void merge_and_clean(std::vector<Particle>& p, double* s)
{
    auto s_tmp = sfc(p);
    array_merge(s,s_tmp,GridNum);
    std::vector<Particle>().swap(p);
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

void snap_push_back(std::string ifname, std::vector<Particle>& p)
{
    std::ifstream ifs(ifname, std::ios_base::binary);
    if(!ifs){
        std::cout << "open file " << ifname << "with error, Abort" << std::endl;
        std::terminate();
    }
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
    std::cout << "p.size(): " << p.size() << std::endl;
}
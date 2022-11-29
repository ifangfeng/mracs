#include"mracs.h"

double* grid_cic(std::vector<Particle>& p, double R, int JCIC);

int main()
{
    read_parameter();
    //auto p10 = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto p20 = read_in_Halo_4vector("/data0/BigMDPL/BigMDPL_halo.bin");
    
    std::vector<Particle> p1, p2;  

    //for(size_t i = 0; i < p10.size(); i += 1) p1.push_back({p10[i].x, p10[i].y, p10[i].z, 1.});
    //std::vector<Particle>().swap(p10);
    //std::cout << "dm: " << p1.size() << std::endl;

    const double M_min {2e12};
    for(size_t i = 0; i < p20.size(); ++i) if(p20[i].weight > M_min) p2.push_back({p20[i].x, p20[i].y, p20[i].z, 1.});
    std::vector<Particle>().swap(p20);
    std::cout << "halo: " << p2.size() << std::endl;
    //#################################################

    const int JCIC {10};
    auto rho_cic = grid_cic(p2, Radius, JCIC);
    auto sum = array_sum(rho_cic,1L<<JCIC*3);
    //double exp = p2.size() 
    std::cout << "ave= " << sum / (1L<<JCIC*3) << std::endl;
    for(int i = 0; i < (1<<JCIC*3); i+=1000000) std::cout << rho_cic[i] << ", "; std::cout << std::endl;



}

double* grid_cic(std::vector<Particle>& p, double R, int JCIC)
{
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    const int64_t Len {1L<<JCIC};
    const int64_t Num {1L<<JCIC*3};
    auto a = new double[Num];
    for(size_t i = 0; i < Num; ++i) a[i] = 0;

    const double Rgrid = R / SimBoxL * Len;
    const double incrmt = 3./(4 * M_PI * pow(R, 3));
    std::cout << "incrmt= " << incrmt << std::endl;
    const int lmin = -floor(Rgrid);
    const int lmax = floor(Rgrid) + 1;
    const int lbox = lmax*2;
    int64_t Numtemp = pow(lbox, 3);
    
    auto temparr = new double[Numtemp];
    for(size_t i = 0; i < Numtemp; ++i) temparr[i] = 0;
    for(int i = lmin; i <= lmax; ++i)
        for(int j = lmin; j <= lmax; ++j)
            for(int k = lmin; k <= lmax; ++k)
                if(i * i + j * j + k * k  < Rgrid * Rgrid)
                    temparr[(i - lmin) * lbox * lbox + (j - lmin) * lbox + (k - lmin)] = incrmt;

    for(size_t n = 0; n < p.size(); ++n)
    {
        int xc = floor(p[n].x / SimBoxL * Len + 0.5);
        int yc = floor(p[n].y / SimBoxL * Len + 0.5);
        int zc = floor(p[n].z / SimBoxL * Len + 0.5);
        for(int i = lmin; i <= lmax; ++i)
            for(int j = lmin; j <= lmax; ++j)
                for(int k = lmin; k <= lmax; ++k)
                    a[((xc + i) & (Len - 1)) * Len * Len + ((yc + j) & (Len - 1)) * Len + (zc + k) & (Len - 1)] 
                    += temparr[(i - lmin) * lbox * lbox + (j - lmin) * lbox + (k - lmin)];
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 6 cic grid = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    delete[] temparr;
    return a;
}
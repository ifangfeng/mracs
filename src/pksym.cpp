#include"mracs.h"

void print_and_clear(double* Pk_grid, double* Pk_gridFree, double* s, int a, int n);
void density_power(int a, int n, std::vector<Particle> p);

int main()
{
    read_parameter();
    //auto p  = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    //auto p0 = default_random_particle(SimBoxL,p.size());

    
    //density_power(1,4,p);
    //density_power(1,10,p);
    density_power(0,5,p);
    //density_power(1,4,p);
    //density_power(0,2,p);
    //density_power(0,3,p);
    //density_power(0,5,p);

}

void print_and_clear(double* Pk_grid, double* Pk_gridFree, double* s, int a, int n)
{
    std::cout << "=====================a: " << a << " , n: " << n << "==========================" << std::endl;
    //for(int i = 0; i < GridLen/2+1; ++i) std::cout << TWOPI*i/SimBoxL             << ", "; std::cout << '\n' ;
    for(int i = 0; i < GridLen/2 +1; ++i) std::cout << Pk_grid[i]            << ", "; std::cout << '\n';
    for(int i = 0; i < GridLen/2 +1; ++i) std::cout << Pk_gridFree[i]        << ", "; std::cout << '\n';
    //for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_grid[i]*i*i*i      << ", "; std::cout << '\n' << '\n';
    //for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_gridFree[i]*i*i*i  << ", "; std::cout << '\n';
    delete[] Pk_grid;
    delete[] Pk_gridFree;
    delete[] s;
}

void density_power(int a, int n, std::vector<Particle> p)
{
    force_base_type(a, n);
    auto s = sfc(p);

    //force_kernel_type(1);
    //auto w = wfc(2,0);
    //#pragma omp parallel for
    //for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    //{
    //    sc[i][0] *= w[i];
    //    sc[i][1] *= w[i]; 
    //}

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto sc= sfc_r2c(s,false);

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference FFT3d    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    double* Pk_grid     = densityPowerFFT(sc,false);
    double* Pk_gridFree = densityPowerDWT(sc,false);
    print_and_clear(Pk_grid, Pk_gridFree, s, a, n);
}
#include"mracs.h"

void print_and_clear(double* Pk_grid, double* Pk_gridFree, double* s, int a, int n);
void density_power(int a, int n, std::vector<Particle> p);

int main()
{
    read_parameter();
    auto p = read_in_Halo_3vector(DataDirec);

    density_power(0,1,p);
    density_power(0,5,p);
    density_power(1,4,p);
    density_power(1,10,p);

}

void print_and_clear(double* Pk_grid, double* Pk_gridFree, double* s, int a, int n)
{
    std::cout << "=====================a: " << a << " , n: " << n << "==========================" << std::endl;
    for(int i = 0; i < GridLen/2+1; ++i) std::cout << TWOPI*i/SimBoxL             << ", "; std::cout << '\n' << '\n';
    for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_grid[i]            << ", "; std::cout << '\n' << '\n';
    for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_gridFree[i]        << ", "; std::cout << '\n' << '\n';
    for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_grid[i]*i*i*i      << ", "; std::cout << '\n' << '\n';
    for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_gridFree[i]*i*i*i  << ", "; std::cout << '\n';
    delete[] Pk_grid;
    delete[] Pk_gridFree;
    delete[] s;
}

void density_power(int a, int n, std::vector<Particle> p)
{
    force_base_type(a, n);
    auto s = sfc(p);
    auto sc= sfc_r2c(s,false);
    double* Pk_grid     = densityPowerFFT(sc);
    double* Pk_gridFree = densityPowerDWT(sc);
    print_and_clear(Pk_grid, Pk_gridFree, s, a, n);
}
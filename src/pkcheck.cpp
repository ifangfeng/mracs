#include"mracs.h"

int main()
{
    read_parameter();
    force_kernel_type(1);
    std::vector<double> R_vec{8,10,20,30,40,50};
    std::vector<double> var, var_th;
    auto p = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto s = sfc(p);
    auto sc= sfc_r2c(s);
    for(auto R : R_vec){
        auto w = wfc(R,0);
        auto c = convol_c2r(sc,w);
        var.push_back(inner_product(c,c,GridNum)/pow(p.size()*4./3*M_PI*pow(R/SimBoxL,3),2)/GridNum-1);
    }

    auto Pk = densityPowerDWT(s);
    double itg{0};
    const double deltak{TWOPI/GridLen};
    double fr[GridLen/2];
    for(auto R : R_vec){
        for(int i = 0; i < GridLen/2; ++i) fr[i] = 3/(pow(i*R*deltak,3))*(sin(i*R*deltak)-i*R*deltak*cos(i*R*deltak));fr[0]=1;
        for(int i = 0; i < GridLen/2; ++i)
        {
            
            itg += i * i * fr[i] * fr[i] * Pk[i];
        }
        var_th.push_back(itg);
    }
    std::cout << "R: " << std::endl; for(auto x : R_vec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "sigma: " << std::endl; for(auto x : var) std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "sigma_th : " << std::endl; for(auto x : var_th) std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "sigma/Sigma_th: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i])/sqrt(var_th[i]) << ", "; std::cout << std::endl;
}
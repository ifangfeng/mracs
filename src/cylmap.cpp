// cylinder map
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.rsd.bin");
    auto p_non = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    force_kernel_type(4);

    const int NumH{20}, NumD{20}; // shape: H/D; 
    const double Hmin{10}, Dmin{10}, Hmax{100}, Dmax{100};
    std::vector<double> Hvec, Dvec, var, var_non, xvec, yvec;
    for(int i = 0; i <= NumH; ++i) Hvec.push_back(Hmin + i*(Hmax-Hmin)/NumH);
    for(int i = 0; i <= NumD; ++i) Dvec.push_back(Dmin + i*(Dmax-Dmin)/NumD);

    auto s = sfc(p); std::cout << "raw var: " << inner_product(s,s,GridVol) * GridVol/pow(p.size(), 2) - 1 << std::endl;
    auto sc= sfc_r2c(s);
    auto s_non = sfc(p_non); std::cout << "raw var: " << inner_product(s_non,s_non,GridVol) * GridVol/pow(p_non.size(), 2) - 1 << std::endl;
    auto sc_non= sfc_r2c(s_non);
    for(auto H : Hvec)
        for(auto D : Dvec)
        {
            auto w = wfc(D/2,H);
            auto c = convol_c2r(sc,w);
            auto c_non = convol_c2r(sc_non,w);
            var.push_back(inner_product(c,c,GridVol) * GridVol/pow(p.size(), 2) - 1);
            var_non.push_back(inner_product(c_non,c_non,GridVol) * GridVol/pow(p_non.size(), 2) - 1);
            xvec.push_back(D); yvec.push_back(H);
            delete[] w;
            delete[] c;
            delete[] c_non;
        }
    std::cout << "H: "   << std::endl;for(auto x : yvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "D: "   << std::endl;for(auto x : xvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "Sigma: " << std::endl;for(auto x : var)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "Sigma_non: " << std::endl;for(auto x : var_non)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "sigma-sigma_non: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i]) - sqrt(var_non[i]) << ", "; std::cout << std::endl;
    std::cout << "sigma/sigma_non: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i]) / sqrt(var_non[i]) << ", "; std::cout << std::endl;

}
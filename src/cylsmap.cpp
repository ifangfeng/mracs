// cylinder shape map
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.rsd.bin");
    auto p_non = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    force_kernel_type(4);

    const int NumShape{20}, NumScale{20}; // shape: H/D; scale: (H*D^2)^(1/3)
    const double SPmin{0.05}, SPmax{20}, SLmin{10}, SLmax{100}; // Volume = Pi/4*SL^3,
    std::vector<double> SPvec, SLvec, var, var_non, xvec, yvec, Dvec, Hvec;
    for(int i = 0; i <= NumShape; ++i) SPvec.push_back(SPmin * pow(SPmax/SPmin,static_cast<double>(i)/NumShape));
    for(int i = 0; i <= NumScale; ++i) SLvec.push_back(SLmin * pow(SLmax/SLmin,static_cast<double>(i)/NumScale));

    auto s = sfc(p); std::cout << "raw var: " << inner_product(s,s,GridVol) * GridVol/pow(p.size(), 2) - 1 << std::endl;
    auto sc= sfc_r2c(s);
    auto s_non = sfc(p_non); std::cout << "raw var: " << inner_product(s_non,s_non,GridVol) * GridVol/pow(p_non.size(), 2) - 1 << std::endl;
    auto sc_non= sfc_r2c(s_non);
    for(auto SP : SPvec)
        for(auto SL : SLvec)
        {
            double D = SL * pow(SP,-1./3);
            double H = SL * pow(SP, 2./3);
            auto w = wfc(D/2,H);
            auto c = convol_c2r(sc,w);
            auto c_non = convol_c2r(sc_non,w);
            var.push_back(inner_product(c,c,GridVol) * GridVol/pow(p.size(), 2) - 1);
            var_non.push_back(inner_product(c_non,c_non,GridVol) * GridVol/pow(p_non.size(), 2) - 1);
            xvec.push_back(SL); yvec.push_back(SP); Hvec.push_back(H); Dvec.push_back(D);
            delete[] w;
            delete[] c;
            delete[] c_non;
        }
    std::cout << "H: "   << std::endl;for(auto x : Hvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "D: "   << std::endl;for(auto x : Dvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "Scale: "   << std::endl;for(auto x : xvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "Shape: "   << std::endl;for(auto x : yvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "Sigma: " << std::endl;for(auto x : var)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "Sigma_non: " << std::endl;for(auto x : var_non)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "sigma-sigma_non: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i]) - sqrt(var_non[i]) << ", "; std::cout << std::endl;
    std::cout << "sigma/sigma_non: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i]) / sqrt(var_non[i]) << ", "; std::cout << std::endl;


}
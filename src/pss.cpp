// power spectrum to second order statistics
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto s = sfc(p);
    auto sc= sfc_r2c(s,false);
    auto Pk = densityPowerDWT(sc);

    const int NumH{3}, NumD{1}; // shape: H/D; 
    const double Hmin{10}, Dmin{100}, Hmax{100}, Dmax{100};
    std::vector<double> Hvec, Dvec, var, var_th, var_th2, xvec, yvec;
    for(int i = 0; i <= NumH; ++i) Hvec.push_back(Hmin + i*(Hmax-Hmin)/NumH);
    for(int i = 0; i < NumD; ++i) Dvec.push_back(Dmin + i*(Dmax-Dmin)/NumD);

    const double deltaK {TWOPI/SimBoxL};
    const int FineConst {10}; // interpolate parameter
    double* fr = new double[GridLen/2 * GridLen/2];
    double fz[GridLen/2];
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    for(auto D : Dvec)
        for(auto H : Hvec)
        {
            for(int i = 0; i < GridLen/2; ++i)
            {
                fz[i] = sin(i*H*deltaK/2)/(i*H*deltaK/2);fz[0]=1;
            }
            for(int i = 0; i < GridLen/2; ++i)
                for(int j = 0; j < GridLen/2; ++j)
                {
                    fr[i*GridLen/2 + j] = std::cyl_bessel_j(1,sqrt(i*i +j*j)*D/2*deltaK)/(sqrt(i*i +j*j)*D/2*deltaK/2);
                }fr[0]=1;
            double itg{0};
            for(int i = 0; i < GridLen/3; ++i)
                for(int j = 0; j < GridLen/3; ++j)
                    for(int k = 0; k < GridLen/3; ++k)
                {
                    
                        itg += fr[i*GridLen/2 + j] * fr[i*GridLen/2 + j] * fz[k] * fz[k] * Pk[i*GridLen/2*GridLen/2 + j*GridLen/2 +k];
                    
                }
            //itg *= pow(deltaK,1);
            var_th.push_back(itg);
        }
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    for(auto D : Dvec)
        for(auto H : Hvec)
        {
            auto w = wfc(D/2,H);
            auto c = convol_c2r(sc,w);
            var.push_back(inner_product(c,c,GridVol) * GridVol/pow(p.size(), 2) - 1);
            auto pk = densityPowerDWT(sc);
            double itg{0};
            for(size_t i = 0; i < pow((GridLen/2),3); ++i) itg += pk[i];
            var_th2.push_back(itg);
            xvec.push_back(D); yvec.push_back(H);
            delete[] w;
            delete[] c;
            delete[] pk;
        }
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 1 th    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;
    std::cout << "Time difference 2 DWt    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count()
    << "[ms]" << std::endl;

    std::cout << "H: "   << std::endl;for(auto x : yvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "D: "   << std::endl;for(auto x : xvec) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "Sigma: " << std::endl;for(auto x : var)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "Sigma_th: " << std::endl;for(auto x : var_th)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "Sigma_th2: " << std::endl;for(auto x : var_th2)  std::cout << sqrt(x) << ", "; std::cout << std::endl;
    std::cout << "sigma/Sigma_th: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i])/sqrt(var_th[i]) << ", "; std::cout << std::endl;
    std::cout << "sigma/Sigma_th2: " << std::endl;for(int i = 0; i < var.size(); ++i)  std::cout << sqrt(var[i])/sqrt(var_th2[i]) << ", "; std::cout << std::endl;

    //std::cout << "=====================a: " << BaseType << " , n: " << phiGenus << "==========================" << std::endl;
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << TWOPI*i/SimBoxL             << ", "; std::cout << '\n' << '\n';
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << Pk[i]            << ", "; std::cout << '\n' << '\n';
    //for(int i = 0; i < GridLen/2 + 1; ++i) std::cout << Pk[i]*i*i*i      << ", "; std::cout << '\n' << '\n';


}

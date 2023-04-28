// variance and corelation with cylinder filter
#include"mracs.h"

int main(){
    read_parameter();
    
    auto p = read_in_Halo_3vector(DataDirec);
    auto s = sfc(p);
    force_kernel_type(4);
    double H{10};
    auto w = wft(Radius,H);
    auto s1= convol3d(s,w);
    force_kernel_type(0);
    std::vector<double> v_xi, v_theta, v_radii;
    const double R0{5}, R1{50}, N_RADII{5}, N_THETA{5};
    for(int i = 0; i < N_THETA; ++i) v_theta.push_back(M_PI/2 * static_cast<double>(i+1)/N_THETA);
    for(int i = 0; i < N_RADII; ++i) v_radii.push_back(R0 * pow((R1/R0), static_cast<double>(i)/N_RADII));
        
    auto s = sfc(p);
    auto sc = sfc_r2c(s);
    force_kernel_type(3);

    for(int i = 0; i < N_RADII; ++i)
        for(int j = 0; j < N_THETA; ++j)
        {
            auto w = wfc(v_radii[i], v_theta[j]);
            auto c = convol_c2r(sc, w);
            v_xi.push_back(inner_product(s, c, GridVol)*GridVol/pow(p.size(), 2) - 1);
            delete[] c;
            delete[] w;
        }

    std::cout << "theta/Pi: ";
    for(auto x : v_theta) std::cout << x/M_PI << ", ";
    std::cout << std::endl;
    for(int i = 0; i < N_RADII; ++i)
    {   
        std::cout << "R=" << v_radii[i] << ": ";
        for(int j =0; j < N_THETA; ++j)
            std::cout << v_xi[i*N_THETA + j] << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    //auto w1= wft()
    //auto Pk = densityPowerDWT(s);
}
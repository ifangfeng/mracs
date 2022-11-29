#include"mracs.h"

#define R0 0.1           // Mpc/h
#define R1 150.1           // Mpc/h
#define NUMTEST 100

int main()
{
    read_parameter();
    std::vector<Particle> p = read_in_DM_3vector(DataDirec);
    /*for(int i = 0; i < 10; ++i)
    {
        std::cout << p[i].x << ", " << p[i].y << ", " << p[i].z << ", " << p[i].weight << std::endl;
    }*/
    int count {0};
    for(int i = 0; i < p.size(); ++i)
    {
        if(p[i].x > 1000 || p[i].x < 0 || p[i].y > 1000 || p[i].y < 0 || p[i].z > 1000 || p[i].z < 0 || p[i].weight != 1.)
        //std::cout << i << ": "<< p[i].x << ", " << p[i].y << ", " << p[i].z << ", " << p[i].weight << std::endl;
        count++;
    }
    std::cout << count << std::endl;
    /*
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);
    std::vector<Particle> p0; for(size_t i = 0; i < p.size(); ++i) p0.push_back({u(e), u(e), u(e), p[i].weight});
    
    auto s   = sfc(p);
    auto sc  = sfc_r2c(s);
    auto s0  = sfc(p0);
    auto sc0 = sfc_r2c(s0);

    std::vector<double> r_log, xi_r, xi_feng, xi_random, xi_r_H, xi_r_LS, xi_r_ph, xi_r_DP;
    for(int i = 0; i < NUMTEST; ++i){
        std::cout << "|" << std::setw(3) << i << "| th point: ";
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        auto c0 = convol_c2r(sc0, w);
        double dd  = inner_product(s, c, GridNum);
        double dr1 = inner_product(s0,c, GridNum);
        double dr2 = inner_product(s,c0, GridNum);
        double rr  = inner_product(s0,c0,GridNum);
        double dr  = (dr1 + dr2)/2;
        xi_r_DP.push_back(dd/dr-1);
        xi_r_ph.push_back(dd/rr-1);
        xi_r_H.push_back(dd*rr/dr/dr-1);
        xi_r_LS.push_back((dd-2*dr+rr)/rr);
        xi_r.push_back(dd * GridNum/pow(p.size(), 2) - 1);
        xi_feng.push_back((dd - rr ) * GridNum/pow(p.size(), 2));
        xi_random.push_back(rr * GridNum/pow(p.size(), 2) - 1);
        delete[] w;
        delete[] c;
        delete[] c0;
    }

    for(auto x : r_log)  std::cout << x << ", "; std::cout << "\nxi_r:\n";
    for(auto x : xi_r)   std::cout << x << ", "; std::cout << std::endl;
    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_r[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_r_LS:\n";
//
//    for(auto x : xi_r_DP)  std::cout << x << ", "; std::cout << std::endl;
//    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_r_DP[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_r_H:\n";
//
//    for(auto x : xi_r_H)  std::cout << x << ", "; std::cout << std::endl;
//    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_r_H[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_r_ph:\n";
//
//    for(auto x : xi_r_ph)  std::cout << x << ", "; std::cout << std::endl;
//    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_r_ph[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_r_LS:\n";
//
    for(auto x : xi_r_LS)  std::cout << x << ", "; std::cout << std::endl;
    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_r_LS[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_feng:\n";

    for(auto x : xi_feng)  std::cout << x << ", "; std::cout << std::endl;
    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_feng[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_random:\n";

    for(auto x : xi_random)  std::cout << x << ", "; std::cout << std::endl;
    for(int i = 0; i < NUMTEST; ++i) std::cout << xi_random[i] * r_log[i] * r_log[i] << ", "; std::cout << "\nxi_random:\n";
    */
}
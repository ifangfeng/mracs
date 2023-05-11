//compare different Correlation result
#include"mracs.h"

int main(){
    read_parameter();
    std::string PREFIX {"output/xi_cp_J" + std::to_string(Resolution) + "_"};
    std::vector<std::string> vofn {"PH","DP","HM","LS","DF1","DF2","DF3","random"};
    std::vector<std::ofstream> vofs;
    for(auto ofn : vofn) 
        vofs.push_back(static_cast<std::ofstream>(PREFIX + ofn));
    const double R0{0.1}, R1{150};
    const int NUMTEST{100};
    std::vector<double> r_log;
    std::ofstream ofs{static_cast<std::ofstream>(PREFIX + "rbin")};
    std::vector<std::vector<double>> vxi;
    vxi.resize(vofn.size());

    auto p  = read_in_DM_3vector(DataDirec);
    auto p0 = default_random_particle(SimBoxL,p.size());

    auto s  = sfc(p);
    auto s0 = sfc(p0);
    auto sc = sfc_r2c(s);
    auto sc0= sfc_r2c(s0);

    auto us = new double[GridVol];
    for(size_t i = 0; i < GridVol; ++i) us[i] = s[i] - s0[i];
    auto usc= sfc_r2c(us);
    
    force_kernel_type(0);
    for(int i = 0; i < NUMTEST; ++i){
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        auto c0 = convol_c2r(sc0, w);
        auto uc = convol_c2r(usc, w);
        double dt  = inner_product(us,uc,GridVol);
        double dd  = inner_product(s, c, GridVol);
        double dr1 = inner_product(s0,c, GridVol);
        double dr2 = inner_product(s,c0, GridVol);
        double rr  = inner_product(s0,c0,GridVol);
        double dr  = (dr1 + dr2)/2;

        vxi[0].push_back({dd/rr-1});                                // PH
        vxi[1].push_back({dd/dr-1});                                // DP
        vxi[2].push_back({dd*rr/(dr*dr)-1});                        // HM
        vxi[3].push_back({(dd-2*dr+rr)/rr});                        // LS
        vxi[4].push_back({dt * GridVol/pow(p.size(), 2)});          // DF1:DT=(d-r)^2/<rr>
        vxi[5].push_back({dd * GridVol/pow(p.size(), 2) - 1});      // DF2:DF=dd/<rr> - 1
        vxi[6].push_back({(dd - rr) * GridVol/pow(p.size(), 2)});   // DF3:FF=(dd-rr)/<rr>
        vxi[7].push_back({rr * GridVol/pow(p.size(), 2) - 1});      // random
        delete[] w;
        delete[] c;
        delete[] c0;
        delete[] uc;
    }

    for(auto r : r_log) ofs << r << ", ";
    for(int i = 0; i < vxi.size(); ++i){
        auto x = std::move(vofs[i]);
        for(int j = 0; j < r_log.size(); ++j) x << vxi[i][j] << ", ";
    }
}

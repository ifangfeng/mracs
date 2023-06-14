// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    std::vector<double> we{-8.81084, 0.0134186, 0.751282, 1}; //weight
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};

    std::ifstream ifs {ifname};
    std::vector<int> envi;int temp{0}; char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) std::cout << "halo size matched! continue\n"; else std::terminate();
    std::vector<double> npartenvi(4,0);
    for(auto x : envi) npartenvi[x]++;
    double sum{0};
    for(auto x : npartenvi) {std::cout << x << ", ";sum+=x;} std::cout << sum << std::endl;
    std::vector<Particle> hlw;
    for(size_t i = 0; i < envi.size(); ++i){
        hlw.push_back({hl[i].x,hl[i].y,hl[i].z,we[envi[i]]*npartenvi[envi[i]]});
    }

    const double R0{1},R1{125};
    const int NUMTEST{20};
    std::vector<double> r_log; for(int i = 0; i < NUMTEST; ++i) r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
    std::vector<double> cross;
    
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc = sfc_r2c(sfc(hlw),true);

    force_kernel_type(1);
    for(auto r : r_log){
        auto w = wfc(r,0);
        auto c1 = convol_c2r(sc,w);
        auto c0 = convol_c2r(sc_dm,w);
        double hm = inner_product(c0,c1,GridVol)*GridVol/dm.size()/hlw.size() -1;
        double hh = inner_product(c1,c1,GridVol)*GridVol/hlw.size()/hlw.size() -1;
        double mm = inner_product(c0,c0,GridVol)*GridVol/dm.size()/dm.size() -1;
        cross.push_back(hm/sqrt(hh*mm));
        delete[] c0;
        delete[] c1;
        delete[] w;
    }
    std::ofstream ofs {"output/ccc_in_real_space_PCA_GSR3.dat"};
    for(auto x : cross){
        ofs << x << " ";
    }ofs << "\n";
    for(auto r : r_log) ofs << r << " ";
}


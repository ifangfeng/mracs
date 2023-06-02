// *******************************************************
// halo environment split and output as separate files
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    /*
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    double GSR {5};
    force_base_type(0,1);
    force_kernel_type(2);
    auto sc = sfc_r2c(sfc(p1),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p2);
    delete[] w;
    delete[] sc;
    */
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR5_halo_Mcut2e12.txt"};
    std::ifstream ifs {ifname};
    std::vector<int> envi;int temp{0}; char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) std::cout << "halo size matched! continue\n"; else std::terminate();

    std::vector<Particle> vd,st,fl,kt;

    for(int i = 0; i < envi.size(); ++i){
        if (envi[i] == 0) vd.push_back(hl[i]);
        else if (envi[i] == 1) st.push_back(hl[i]);
        else if (envi[i] == 2) fl.push_back(hl[i]);
        else if (envi[i] == 3) kt.push_back(hl[i]);
    }

    std::ofstream ofs_vd {"output/envi_vd_GSR5_halo_Mcut2e12.txt"};
    for(auto x : vd) ofs_vd << x << ", ";
}


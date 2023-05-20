// *******************************************************
// halo environment split and sample for scatter plot
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    std::string GSR {"_GSR5"};
    std::string ifname {"output/envi_J10" + GSR + "_halo_Mcut2e12.txt"};
    std::ifstream ifs {ifname};

    std::vector<char> envi; char temp{0}, comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    auto p0 = generate_random_particle(20,SimBoxL,0);
    std::vector<double> vecR {50};

    auto s1 = sfc(dm);
    auto s2 = sfc(hl);
    auto sc1 = sfc_r2c(s1);
    auto sc2 = sfc_r2c(s2);


    std::string prefix {"output/envi_split_scatter_algo2_"};
    std::string suffix {"_halo_Mcut2e12.txt"};
    std::vector<double> result, nvd, nst, nfl, nkt;
    
    force_kernel_type(1);
    for(auto R : vecR){
        auto w = wfc(R,0);
        auto c1 = convol_c2r(sc1,w);
        auto c2 = convol_c2r(sc2,w);
        auto n1 = project_value(c1,p0);
        auto n2 = project_value(c2,p0);
        
        double cicexpect1 = static_cast<double>(dm.size())/GridVol;
        std::string ofname_dm = prefix + "dm_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_dm {ofname_dm};
        for(size_t i = 0; i < p0.size(); ++i) ofs_dm << n1[i]/cicexpect1 << ", ";

        double cicexpect2 = static_cast<double>(hl.size())/GridVol;
        for(size_t i = 0; i < p0.size(); ++i) result.push_back(n2[i]/cicexpect2);
        std::string ofname_hl = prefix + "hl_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_hl {ofname_hl};
        for(auto x : result) ofs_hl << x << ", ";

        for(size_t i = 0; i < p0.size(); ++i)
        {
            if(envi[i] == '0') nvd.push_back(result[i]);
            else if(envi[i] == '1') nst.push_back(result[i]);
            else if(envi[i] == '2') nfl.push_back(result[i]);
            else if(envi[i] == '3') nkt.push_back(result[i]);
        }

        std::string ofname_vd = prefix + "vd_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::string ofname_vdx = prefix + "vd_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_vd {ofname_vd};
        std::ofstream ofs_vdx {ofname_vd};
        for(auto x : nvd) ofs_vd << x << ", ";
        std::string ofname_st = prefix + "st_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_st {ofname_st};
        for(auto x : nst) ofs_st << x << ", ";
        std::string ofname_fl = prefix + "fl_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_fl {ofname_fl};
        for(auto x : nfl) ofs_fl << x << ", ";
        std::string ofname_kt = prefix + "kt_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
        std::ofstream ofs_kt {ofname_kt};
        for(auto x : nkt) ofs_kt << x << ", ";

        std::vector<double>().swap(result);
        std::vector<double>().swap(nvd);
        std::vector<double>().swap(nst);
        std::vector<double>().swap(nfl);
        std::vector<double>().swap(nkt);
    }
}

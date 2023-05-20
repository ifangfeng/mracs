// *******************************************************
// halo environment split and sample for scatter plot
// *******************************************************
#include"mracs.h"

std::string GSR {"_GSR5"};

int main(){
    read_parameter();
    std::string ifname {"output/envi_J10" + GSR + "_halo_Mcut2e12.txt"};
    std::string ifnameGrid {"output/envi_J10" + GSR + "_halo_Mcut2e12_grid.txt"};
    std::ifstream ifs {ifname}, ifsg {ifnameGrid};

    std::vector<char> envi, envig; char temp{0}, comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    while(ifsg >> temp >> comma) envig.push_back(temp);

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    if(envi.size() == hl.size()) std::cout << "size matched, continue" << std::endl;
    else {std::cout << "loading data with error, Abort\n"; return 0;}
    if(envig.size() == 1024*1024*1024) std::cout << "size matched, continue" << std::endl;
    else {std::cout << "loading data with error, Abort\n"; return 0;}

    std::vector<Particle> voids,sheets,filaments,knots;
    for(size_t i = 0; i < hl.size(); ++i){
        if(envi[i] == '0') voids.push_back(hl[i]);
        else if(envi[i] == '1') sheets.push_back(hl[i]);
        else if(envi[i] == '2') filaments.push_back(hl[i]);
        else if(envi[i] == '3') knots.push_back(hl[i]);
    }
    std::vector<int64_t> voidsg,sheetsg,filamentsg,knotsg;
    for(int64_t i = 0; i < envig.size(); ++i){
        if(envig[i] == '0') voidsg.push_back(i);
        else if(envig[i] == '1') sheetsg.push_back(i);
        else if(envig[i] == '2') filamentsg.push_back(i);
        else if(envig[i] == '3') knotsg.push_back(i);
    }

    auto p0 = generate_random_particle(20,SimBoxL,0);
    std::vector<double> vecR {30};

    std::vector<std::vector<Particle>> vecP {voids, sheets, filaments, knots, hl, dm};
    std::vector<std::vector<int64_t>> vecPg {voidsg, sheetsg, filamentsg, knotsg};
    std::vector<double*> svec, svecg;
    for(auto x : vecP) svec.push_back(sfc(x));
    for(auto x : vecPg) svecg.push_back(sfc_grid_coordinate(x));

    std::string prefix {"output/envi_split_scatter_algo3_"};
    std::vector<std::string> fname {"vd","st","fl","kt","hl","dm"};
    std::string suffix {"_halo_Mcut2e12.txt"};
    
    force_kernel_type(1);
    for(auto R : vecR){
        std::cout << "processing R" << R << ":\n";
        auto w = wfc(R,0);
        for(int i = 0; i < vecP.size(); ++i){
            auto c = convol3d(svec[i],w);
            auto n = project_value(c,p0);
            std::string ofname = prefix + fname[i] +"_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
            std::ofstream ofs {ofname};
            if(i > 3){   
                double cicexpect = static_cast<double>(vecP[i].size())/GridVol;
                for(size_t j = 0; j < p0.size(); ++j) ofs << n[j]/cicexpect << ", ";
            }
            else{
                auto cg = convol3d(svecg[i],w);
                auto ng = project_value(cg,p0);
                for(int64_t j = 0; j < p0.size(); ++j) ofs << (n[j]*vecPg[i].size())/(ng[j]*vecP[i].size()) << ", ";
                delete[] cg;
                delete[] ng;
            }
            delete[] c;
            delete[] n;
        }
    }
}

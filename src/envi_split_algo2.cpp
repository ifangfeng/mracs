// *******************************************************
// halo environment split and sample for scatter plot
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    auto p0 = generate_random_particle(20,SimBoxL,0);
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<double> vecR {20,30,40,50,80};

    std::vector<std::string> ofd {"_x","_y"};
    std::vector<std::string> ofns {"vd","st","fl","kt"};

    std::string prefix {"output/envi_split_scatter_algo2_"};
    std::string suffix {"_GSR8_halo_Mcut2e12.txt"};


    auto envi = environment(dm,8,p0);

    force_base_type(0,5);
    force_kernel_type(1);
    auto sc1 = sfc_r2c(sfc(dm),true);
    auto sc2 = sfc_r2c(sfc(hl),true);

    for(auto R : vecR){
        auto w = wfc(R,0);
        auto n1 = project_value(convol_c2r(sc1,w), p0, true);
        auto n2 = project_value(convol_c2r(sc2,w), p0, true);

        std::string RJ = "_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution);

        double* nt = nullptr;
        double cic {0};
        for(auto dim : ofd){
            if(dim == "_x") {
                nt = n1; 
                cic = static_cast<double>(dm.size())/GridVol;} 
            else if(dim == "_y") {
                nt = n2; 
                cic = static_cast<double>(hl.size())/GridVol;}
            for(auto name : ofns){
                int marker{0};
                std::string ofname = prefix + name + dim + RJ + suffix;
                std::ofstream ofs {ofname};
                if(name == "st") marker = 1 ;else if(name == "fl") marker = 2;else if(name == "kt") marker = 3; 
                for(size_t i = 0; i < p0.size(); ++i){
                    if(envi[i] == marker){
                        ofs << nt[i]/cic << ", ";
                    }
                }
            }
            std::string ofname = prefix + ((dim == "_x") ? "dm" : "hl") + RJ + suffix;
            std::ofstream ofs {ofname};
            for(size_t i = 0; i < p0.size(); ++i){
                ofs << nt[i]/cic << ", ";
            }
        }
    }
}

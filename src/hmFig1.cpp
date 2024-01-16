#include"mracs.h"

int main() {
    read_parameter();
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    
    std::vector<int> vec_resol{8,9,10};
    std::vector<std::ofstream*> vec_ofs;
    for(auto x : vec_resol){
        std::ofstream* ofs = new std::ofstream("output/hmFig1_"+RADII+"_J"+std::to_string(x)+".txt");
        vec_ofs.push_back(ofs);
    }

    for(auto ofs : vec_ofs){
        *ofs << "haha\n";
        delete ofs;
    }

    vec_ofs.clear();

    std::cout << std::endl;
}
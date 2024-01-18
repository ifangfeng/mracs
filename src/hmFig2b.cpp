#include"mracs.h"
#include<chrono>

int main() {
    read_parameter();

    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");

    std::vector<int> vec_resol{8,9,10};
    std::vector<std::ofstream*> vec_ofs;
    std::ofstream ofs_cic{"output/hmFig2b.txt"};
    for(auto x : vec_resol){
        std::ofstream* ofs = new std::ofstream("output/hmFig2b_J"+std::to_string(x)+".txt");
        vec_ofs.push_back(ofs);
    }

    std::vector<Particle> p0,p1;
    Offset xx ({450.,450.,500.});
    for(auto x : p){
        if(x.x >= 400 && x.x < 600)
            if(x.y >= 400 && x.y < 600)
                if(x.z > 450 && x.z < 550)
                    p1.push_back(x);
    }
    for(int i = 0; i < 100; ++i)
        for(int j = 0; j < 100; ++j){
            p0.push_back({i+xx.dx,j+xx.dy,xx.dz,1.});
        }
    for(int i = 0; i < vec_ofs.size(); ++i){
        force_resoluton_J(vec_resol[i]);
        auto s = sfc(p);
        auto n = project_value(s,p0,true);
        for(size_t j = 0; j < p0.size(); ++j) *vec_ofs[i] << n[j] << ", ";
        delete[] n;
        delete vec_ofs[i];
    }
    
    for(auto x : p1){
        if(x.x >= 450 && x.x < 550)
            if(x.y >= 450 && x.y < 550)
                if(x.z > 495 && x.z < 505)
                    ofs_cic << x.x -450 << ", " << x.y -450 << ", ";
    }

    vec_ofs.clear();
}
// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto p2 = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    int Nl {1};
    int Nrl {Nl * Nl * Nl}; // number of realization
    int NPW {3};   // number of element in set {<m,h>,<m,m>,<h,h>}
    double Rs {2}; // Gaussian smoothing radius 

    std::vector<std::vector<Particle>> dm(Nrl), hl(Nrl), vd(Nrl), st(Nrl), fl(Nrl), kt(Nrl);
    std::vector<std::vector<std::vector<Particle>>> data {hl,vd,st,fl,kt};
    std::cout << "marker++++++++++++++" << "\n";
        std::cout << "x[0].size():" << data[0][0].size() << ", addr:" << &data[0][0] << "\n";
            std::cout << "hl.size():" << hl[0].size() << ", addr:" << &hl[0] << "\n";
    for(auto x : p1) dm[static_cast<int>(x.x/SimBoxL*Nl) * Nl * Nl + static_cast<int>(x.y/SimBoxL*Nl) * Nl + static_cast<int>(x.z/SimBoxL*Nl)].push_back(x);
    for(auto x : p2) hl[static_cast<int>(x.x/SimBoxL*Nl) * Nl * Nl + static_cast<int>(x.y/SimBoxL*Nl) * Nl + static_cast<int>(x.z/SimBoxL*Nl)].push_back(x);

    for(int i = 0; i < Nrl; ++i){
        auto envi = environment(dm[i],Rs,hl[i]);
        for(int j = 0; j < hl[i].size(); ++j){
            if(envi[j] == 0) vd[i].push_back(hl[i][j]);
            else if(envi[j] == 1) st[i].push_back(hl[i][j]);
            else if(envi[j] == 2) fl[i].push_back(hl[i][j]);
            else if(envi[j] == 3) kt[i].push_back(hl[i][j]);
        }
        if(i == 0){
        std::cout << "dm: " << dm[i].size() << ", hl: " << hl[i].size() << ", vd: " << vd[i].size()
        << ", st: " << st[i].size() << ", fl: " << fl[i].size() << ", kt: " << kt[i].size() << "\n";
        }

        for(int idx = 0;auto x : data){
            if(idx == 0){
            std::cout << x.size() << "\n";
            std::cout << "x[0].size():" << x[0].size() << ", addr:" << &x[0] << "\n";
            std::cout << "hl.size():" << hl[0].size() << ", addr:" << &hl[0] << "\n";
            }
            ++idx;
        }
    }
}
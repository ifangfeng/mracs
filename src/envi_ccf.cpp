// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    int Nl {4};
    int Nrl {Nl * Nl * Nl}; // number of realization
    int NPW {3};   // number of element in set {<m,h>,<m,m>,<h,h>}
    double Rs {5}; // Gaussian smoothing radius 
    
    std::vector<std::string> filename {"hl","vd","st","fl","kt"};
    std::string GSR {"_GSR" + std::to_string((int)Rs)};
    std::string prefix {"output/envi_ccf_"};
    std::string suffix {GSR + "_halo_Mcut2e12.txt"};


    std::vector<std::vector<Particle>> dm(Nrl), hl(Nrl), vd(Nrl), st(Nrl), fl(Nrl), kt(Nrl);
    std::vector<std::vector<std::vector<Particle>>*> data {&hl,&vd,&st,&fl,&kt};
    std::vector<double> column(GridLen,0);
    std::vector<std::vector<double>> pw2d(NPW,column);   // power and cross power of dm and halo
    std::vector<std::vector<std::vector<double>>> pw(data.size(),pw2d);
    std::vector<std::vector<double>> r(data.size(),column);
    std::vector<double *> pt(NPW,nullptr);

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
        auto sc = sfc_r2c(sfc(dm[i]),true);
        pt[0] = densityPowerDWT(sc);
        for(int index = 0; auto x : data){
            auto sc1 = sfc_r2c(sfc((*x)[i]),true);
            pt[1] = densityPowerDWT(sc1);
            pt[2] = crossPowerDWT(sc,sc1);

            for(int j = 0; j < NPW; ++j){
                for(int k = 0; k < GridLen; ++k){
                    pw[index][j][k] += pt[j][k]/Nrl;
                }
            }
            ++index;
        }
    }
    for(int i = 0; i < data.size(); ++i){
        for(int k = 1; k < GridLen; ++k) 
            r[i][k] = pw[i][2][k]/sqrt(pw[i][0][k] * pw[i][1][k]);r[i][0] = 1;
    }

    for(int i = 0; auto x : r){
        std::string ofn = prefix + filename[i] + suffix;
        std::ofstream ofs {ofn};
        for(int k = 0; k < GridLen; ++k){
            ofs << x[k] << ", ";
        }
        ++i;
    }
    std::string ofn = prefix + "kbin" + suffix;
    std::ofstream ofs {ofn};
    for(int i = 0; i < GridLen; ++i)
        ofs << i * TWOPI/(SimBoxL/Nl) << ", "; std::cout << std::endl;

}
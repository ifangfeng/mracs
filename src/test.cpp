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

    std::vector<int64_t> ps;
    for(int64_t i = 0; i < GridVol; ++i)
    {
        if(envig[i] == '1') ps.push_back(i);
    }
    std::cout << static_cast<double>(ps.size())/GridVol << "\n";

    auto s = sfc_grid_coordinate(ps);
    std::cout << array_sum(s,GridVol) << " and " << ps.size() << "\n";
    
}

/*std::vector<int64_t> veci{0,5,8},vecj{5,7,1},veck{8,3,0};
    for(auto i : veci)
        for(auto j : vecj)
            for(auto k : veck)
            {
                int64_t l = i * GridLen * GridLen + j * GridLen + k;
                int64_t ri,rj,rk;
                int64_t bi,bj,bk;
                rk = l&(GridLen-1);
                rj = ((l-rk)&(GridLen*GridLen-1))>>Resolution;
                ri = ((l-rj*GridLen-rk)&(GridLen*GridLen*GridLen-1))>>(Resolution*2);
                bk = l&(GridLen-1);
                bj = (l&((GridLen-1)<<Resolution))>>Resolution;
                bi = (l&((GridLen-1)<<(Resolution*2)))>>(Resolution*2);
                std::cout << i << " " << j << " " << k << " ---> ";
                std::cout << ri << " " << rj << " " << rk << " or ";
                std::cout << bi << " " << bj << " " << bk << '\n';
            }*/
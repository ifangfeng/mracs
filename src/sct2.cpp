// dm and hl scatter
#include"mracs.h"

int main(){
    read_parameter();

    std::vector<double> we{-0.973741, 0.00269408, 0.136522, 0.182164}; //weight
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
        hlw.push_back({hl[i].x,hl[i].y,hl[i].z,we[envi[i]]});//npartenvi[envi[i]]});
    }

    std::vector<std::vector<Particle>*> vec_p{&dm,&hl,&hlw};
    std::vector<double*> vec_n;

    auto w = wfc(Radius,0);
    auto p0 = default_random_particle(SimBoxL,1000*5);
    for(size_t i = 0; i < vec_p.size(); ++i) vec_n.push_back(project_value(convol3d(sfc(*vec_p[i]),w,true), p0, true));

    auto tmp = sfc(hlw);
    auto asum = array_sum(tmp,GridVol);
    std::ofstream ofs {"output/scatter_cp_"+RADII+".txt"};
    for(size_t i = 0; i < p0.size(); ++i) 
    ofs << vec_n[0][i]/dm.size()*GridVol << " " << vec_n[1][i]/hl.size()*GridVol << " " << vec_n[2][i]/asum*GridVol << " ";
    
    
}
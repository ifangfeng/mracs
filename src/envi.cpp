#include"mracs.h"
using namespace std;

int main(){
    read_parameter();
    SimBoxL = 100;
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto p2 = read_in_Halo_4vector("data0/MDPL2/halo_position.bin");
    std::vector<Particle> dm,hl;
    for(auto x : p1)
        if(x.x < 100 && x.y < 100 && x.z < 100)
            dm.push_back({x.x, x.y, x.z, x.weight});
    for(auto x : p2)
        if(x.weight > 2e12 && x.x < 100 && x.y < 100 && x.z < 30)
            hl.push_back({x.x, x.y, x.z, x.weight});
    std::vector<Particle>().swap(p1);
    std::vector<Particle>().swap(p2);
    std::cout << "dm: " << dm.size() << std::endl;
    std::cout << "hl: " << hl.size() << std::endl;

    force_base_type(0,1);
    force_kernel_type(2);
    auto s = sfc(dm);
    auto sc = sfc_r2c(s);
    auto w = wfc(Radius, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,hl);

}
// dm and hl scatter
#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::vector<std::vector<Particle>*> vec_p{&dm,&hl};
    std::vector<double*> vec_n;

    auto w = wfc(Radius,0);
    auto p0 = default_random_particle(SimBoxL,1000*1000);
    for(size_t i = 0; i < vec_p.size(); ++i) vec_n.push_back(project_value(convol3d(sfc(*vec_p[i]),w,true), p0, true));

    std::ofstream ofs {"output/scatter_"+RADII+".txt"};
    for(size_t i = 0; i < p0.size(); ++i) ofs << vec_n[0][i]/dm.size()*GridVol << " " << vec_n[1][i]/hl.size()*GridVol << " ";
    
}
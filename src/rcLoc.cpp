// local density and halo mass Multi-split
#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    // ---access loc dens-----
    auto sc = sfc_r2c(sfc(dm),true);
    auto w  = wfc(Radius,0);
    auto c  = convol_c2r(sc,w);
    auto nprj = project_value(c,hl,true);

    std::vector<double> locDens(hl.size()), vecMass(hl.size());
    for(size_t i = 0; i < hl.size(); ++i) locDens[i] = nprj[i];
    for(size_t i = 0; i < hl.size(); ++i) vecMass[i] = hl[i].weight;

    // ----multi-split with Mass---
    int Mbin{4}, Dbin{4};
    auto nodes_D = nodes_of_proto_sort(locDens,Dbin);
    auto nodes_M = nodes_of_proto_sort(vecMass,Mbin);
    std::vector<std::vector<Particle>*> vpts;
    for(int i = 0; i < Dbin * Mbin; ++i) vpts.push_back(new std::vector<Particle>);
    for(size_t i = 0; i < hl.size(); ++i){
        vpts[classify_index<double>(nodes_D,locDens[i]) * Mbin + classify_index<double>(nodes_M,vecMass[i])]->push_back(hl[i]);
    }
    for(auto x : vpts) print_min_max_and_size(*x);

    force_resoluton_J(7);
    auto sol = optimal_solution(dm,vpts,Radius,false);

    double Npart = hl.size();
    std::cout << "NF:\n";
    for(auto x : vpts) std::cout << x->size() / Npart << ", "; std::cout << std::endl;

}
// covariance of different halo bin 
#include"mracs.h"

double minimum_weight(std::vector<Particle>& v);
double maximum_weight(std::vector<Particle>& v);

int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};


    auto vpts = halo_envi_match_and_split(ifname,hl);

    for(int i = 0; auto x : vpts){
        std::cout << "component " << i << ": min=" << minimum_weight(*x) 
        <<  ", max=" << maximum_weight(*x) << ", size()= " << x->size() << "\n";
        ++i;
    }

    // ----reconstruct and check-------
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_dm = sfc_r2c(sfc(dm),true);

    auto sc_rc = optimal_reconstruct(dm,vpts,Radius,true);
    
    auto wpk = window_Pk(Radius,0);

    auto cc1 = correlation_coefficients(sc_dm,sc_hl,wpk);
    auto cc2 = correlation_coefficients(sc_dm,sc_rc,wpk);
    
    std::cout << "Cross-correlation coefficient:\n";
    std::cout << "[mass picture]   default r_m: " << cc1 << "\n";
    std::cout << "-----------RCST----------r_m: " << cc2 << "\n";
    


    
}


//std::vector<std::vector<Particle>*> envi_and_mass_combine_split()



double minimum_weight(std::vector<Particle>& v)
{
    double min{v[0].weight};
    for(size_t i = 1; i < v.size(); ++i){
        if(v[i].weight < min) min = v[i].weight;
    }

    return min;
}

double maximum_weight(std::vector<Particle>& v)
{
    double max{v[0].weight};
    for(size_t i = 1; i < v.size(); ++i){
        if(v[i].weight > max) max = v[i].weight;
    }

    return max;
}


#include"mracs.h"

int main(){
    read_parameter();

    //auto p = read_in_DM_3vector("../L1380dmposrsd.bin");
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto p0 = default_random_particle(SimBoxL,p.size());
    auto vec_r = linear_scale_generator(20,149.5,50,false);

    auto s = sfc(p); auto sc = sfc_r2c(s,false);
    auto s0= sfc(p0);auto sc0= sfc_r2c(s0,false);

    std::vector<double> corr, corr_DF3, corr_random;

    for(int i = 0; i < vec_r.size(); ++i)
    {
        auto w = wfc(vec_r[i],M_PI/2);
        auto c = convol_c2r(sc,w);
        auto c0= convol_c2r(sc0,w);
        corr.push_back((inner_product(s,c,GridVol)/pow(array_sum(s,GridVol)/GridVol,2)/GridVol -1)*pow(vec_r[i],2));
        corr_DF3.push_back(((inner_product(s,c,GridVol)- inner_product(s0,c0,GridVol)) / pow(array_sum(s,GridVol),2) * GridVol)*pow(vec_r[i],2));
        corr_random.push_back((inner_product(s0,c0,GridVol)/pow(array_sum(s0,GridVol)/GridVol,2)/GridVol -1)*pow(vec_r[i],2));
        delete w;
        delete c;
        delete c0;
    }

    for(auto x : corr) std::cout << x << ", ";std::cout << std::endl;
    for(auto x : corr_DF3) std::cout << x << ", ";std::cout << std::endl;
    for(auto x : corr_random) std::cout << x << ", ";std::cout << std::endl;
    std::cout << "smoothed DF3\n";

    for(int i = 0; i < corr_DF3.size(); ++i){
        if(i == 0){
            std::cout << (corr_DF3[i] + corr_DF3[i+1])/2 << ", ";
        }
        else if(i == corr_DF3.size() - 1)
        {
            std::cout << (corr_DF3[i-1] + corr_DF3[i])/2 << ", ";
        }
        else
            std::cout << (corr_DF3[i-1] + corr_DF3[i] + corr_DF3[i+1])/3 << ", ";
    }std::cout << std::endl;

    
    
}
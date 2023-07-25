#include"mracs.h"

int main(){
    read_parameter();

    std::string ofn {"output/voidrsd.txt"};
    std::ofstream ofs {ofn}; 
    if(!ofs) {std::cout << "opening file with error, abort\n"; std::terminate();}


    auto p = read_in_DM_3vector("../L1380dmposrsd.bin");
    auto vec_r = linear_scale_generator(0.5,149.5,150,false);

    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);
    std::vector<double> corr;

    for(int i = 0; i < vec_r.size(); ++i)
        for(int j = 0; j < vec_r.size(); ++j){
            double r = sqrt(pow(vec_r[i],2) + pow(vec_r[j],2));
            double theta = acos(vec_r[j] / r);
            auto w = wfc(r,theta);
            auto c = convol_c2r(sc,w);
            corr.push_back((inner_product(s,c,GridVol)/pow(array_sum(s,GridVol)/GridVol,2)/GridVol -1)*r*r);
            delete w;
            delete c;
        }

    for(auto x : corr) ofs << x << " ";

    
    
}
// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};

    std::ifstream ifs {ifname};
    std::vector<int> envi;int temp{0}; char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) std::cout << "halo size matched! continue\n"; else std::terminate();
    std::vector<Particle> vd,st,fl,kt;
    for(size_t i = 0; i < envi.size(); ++i){
        if(envi[i] == 0) vd.push_back(hl[i]);
        else if(envi[i] == 1) st.push_back(hl[i]);
        else if(envi[i] == 2) fl.push_back(hl[i]);
        else if(envi[i] == 3) kt.push_back(hl[i]);
    }
    double R{100};
    std::vector<std::vector<Particle>*> data{&vd,&st,&fl,&kt};
    std::vector<int64_t> datasize; for(auto x : data) datasize.push_back((*x).size());
    std::vector<fftw_complex*> vec_sc; for(auto x : data) vec_sc.push_back(sfc_r2c(sfc(*x),true));
    std::vector<double> cr_hl,cr_vd,cr_st,cr_fl,cr_kt;
    std::vector<std::vector<double>*> vec_cross{&cr_hl,&cr_vd,&cr_st,&cr_fl,&cr_kt};
    std::vector<double> covar;

    force_kernel_type(1);
    auto w = wfc(R,0);
    for(int i = 0; i < data.size(); ++i){
        auto c0 = convol_c2r(vec_sc[i],w);
        for(int j = 0; j < data.size(); ++j)
        {
            if(j >= i){
            auto c1 = convol_c2r(vec_sc[j],w);
            covar.push_back(inner_product(c0,c1,GridVol)*GridVol/datasize[j]/datasize[i] -1);
            delete[] c1;
            std::cout << "(i,j)=" << i << ", " << j << std::endl;}
        }
        delete[] c0;
    }
    for(auto x : covar) std::cout << x << ", "; std::cout << std::endl;
}

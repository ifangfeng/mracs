// map of weight in N-dimensional vector space, N is the number of parameter bin
#include"mracs.h"

void halo_envi_match(std::string ifn, std::vector<Particle>& hl, std::vector<Particle>& vd
, std::vector<Particle>& st, std::vector<Particle>& fl, std::vector<Particle>& kt);

int main(){
    read_parameter();

    std::vector<double> vec_R {2,4,8,16,32,64};

    force_resoluton_J(9);
    //===>preprocessing dark matter and halo data
    auto dms = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hls = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    std::vector<Particle> vds,sts,fls,kts;
    halo_envi_match(ifname,hls,vds,sts,fls,kts);
    std::vector<std::vector<Particle>*> vpts {&dms,&hls,&vds,&sts,&fls,&kts};

    auto sc_dms = sfc_r2c(sfc(dms),true);
    auto sc_hls = sfc_r2c(sfc(hls),true);
    auto pk_plus0_s = densityCovarianceArray(sc_dms,sc_hls);
    std::vector<double> ccrs, ccr; // cross-correlation coefficient of different smoothing Radius
    for(auto r : vec_R){
        ccrs.push_back(var_CombinewithKernel(pk_plus0_s,r,0));
    }
    fftw_free(sc_dms);
    fftw_free(sc_hls);
    delete[] pk_plus0_s;

    std::vector<Particle> dm,hl,vd,st,fl,kt;
    std::vector<std::vector<Particle>*> vpt {&dm,&hl,&vd,&st,&fl,&kt};
    SimBoxL /= 4;   // sub box length
    for(int i = 0; i < vpts.size(); ++i){
        for(size_t j = 0; j < (*vpts[i]).size(); ++j){
            if(((*vpts[i])[j].x < SimBoxL) && ((*vpts[i])[j].y < SimBoxL) && ((*vpts[i])[j].z < SimBoxL))
                (*vpt[i]).push_back((*vpts[i])[j]);
        }
        std::vector<Particle>().swap(*vpts[i]);
    }
    for(int i = 0; i < vpts.size(); ++i){
        std::cout << (*vpts[i]).size() << " and sub: " << (*vpt[i]).size() << std::endl;
    }

    force_resoluton_J(7);
    //===>Calculate initial value
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_vd = sfc_r2c(sfc(vd),true);
    auto sc_st = sfc_r2c(sfc(st),true);
    auto sc_fl = sfc_r2c(sfc(fl),true);
    auto sc_kt = sfc_r2c(sfc(kt),true);

    auto pk_plus0 = densityCovarianceArray(sc_dm,sc_hl);
    for(auto r : vec_R){
        ccr.push_back(var_CombinewithKernel(pk_plus0,r,0));
    }
    fftw_free(sc_hl);
    delete[] pk_plus0;

    for(int i = 0; i < vec_R.size(); ++i) std::cout << ccrs[i] << ", ";std::cout << std::endl;
    for(int i = 0; i < vec_R.size(); ++i) std::cout << ccr[i]  << ", ";std::cout << std::endl;


    //===>Here we go!
    const int Nbin{32}; // resoution of parameter space
    double d0,d1,d2,d3; // weight of i_th dimension 

    auto ccc = new double [Nbin * Nbin * Nbin * vec_R.size()];

    for(int i = 0; i < Nbin; ++i)
        for(int j = 0; j < Nbin; ++j)
            for(int k = 0; k < Nbin; ++k){
                d1 = static_cast<double>(i)/Nbin;
                d2 = static_cast<double>(j)/Nbin;
                d3 = static_cast<double>(k)/Nbin;
                d0 = 1 - d1 - d2 - d3;
                for(auto r : vec_R){

                }
            }


}




void halo_envi_match(std::string ifn, std::vector<Particle>& hl, std::vector<Particle>& vd
, std::vector<Particle>& st, std::vector<Particle>& fl, std::vector<Particle>& kt)
{
    std::ifstream ifs {ifn};
    if(!ifs){std::cout << "reading " + ifn + " with error, Abort"; std::terminate();}

    std::vector<int> envi;
    int temp{0};
    char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) 
        std::cout << "halo size matched! continue\n"; 
    else std::terminate();

    for(size_t i = 0; i < envi.size(); ++i){
        if(envi[i] == 0) vd.push_back(hl[i]);
        else if(envi[i] == 1) st.push_back(hl[i]);
        else if(envi[i] == 2) fl.push_back(hl[i]);
        else if(envi[i] == 3) kt.push_back(hl[i]);
    }
    std::vector<int>().swap(envi);
    if(hl.size() != (vd.size() + st.size() + fl.size() + kt.size())) 
        std::cout << "Warning! halo environment subset size not matched\n";
}
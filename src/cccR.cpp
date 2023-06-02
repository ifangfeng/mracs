// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
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

    const double R0{1},R1{125};
    const int NUMTEST{20};
    std::vector<std::vector<Particle>*> data{&hl,&vd,&st,&fl,&kt};
    std::vector<int64_t> datasize; for(auto x : data) datasize.push_back((*x).size());
    std::vector<fftw_complex*> vec_sc; for(auto x : data) vec_sc.push_back(sfc_r2c(sfc(*x),true));
    std::vector<double> r_log; for(int i = 0; i < NUMTEST; ++i) r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
    std::vector<double> cr_hl,cr_vd,cr_st,cr_fl,cr_kt;
    std::vector<std::vector<double>*> vec_cross{&cr_hl,&cr_vd,&cr_st,&cr_fl,&cr_kt};
    
    auto sc_dm = sfc_r2c(sfc(dm),true);

    force_kernel_type(1);
    for(auto r : r_log){
        auto w = wfc(r,0);
        for(size_t i = 0; auto sc : vec_sc){
            auto c1 = convol_c2r(sc,w);
            auto c0 = convol_c2r(sc_dm,w);
            double hm = inner_product(c0,c1,GridVol)*GridVol/dm.size()/datasize[i] -1;
            double hh = inner_product(c1,c1,GridVol)*GridVol/datasize[i]/datasize[i] -1;
            double mm = inner_product(c0,c0,GridVol)*GridVol/dm.size()/dm.size() -1;
            (*vec_cross[i]).push_back(hm/sqrt(hh*mm));
            delete[] c0;
            delete[] c1;
            ++i;
        }
        delete[] w;
    }
    std::ofstream ofs {"output/ccc_in_real_space_GSR3.dat"};
    for(auto x : vec_cross){
        for(int i = 0; i < NUMTEST; ++i) ofs << (*x)[i] << " ";
        ofs << "\n";
    }
    for(auto r : r_log) ofs << r << " ";
}

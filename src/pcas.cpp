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
    double R{100};
    std::vector<std::vector<Particle>*> data{&vd,&st,&fl,&kt};
    std::vector<int64_t> datasize; for(auto x : data) datasize.push_back((*x).size());
    std::vector<fftw_complex*> vec_sc; for(auto x : data) vec_sc.push_back(sfc_r2c(sfc(*x),true));
    std::vector<double> cr_hl,cr_vd,cr_st,cr_fl,cr_kt;
    std::vector<std::vector<double>*> vec_cross{&cr_hl,&cr_vd,&cr_st,&cr_fl,&cr_kt};
    std::vector<double> covar, covar2;

    force_kernel_type(1);
    std::chrono::steady_clock::time_point time0 = std::chrono::steady_clock::now();
    auto w = wfc(R,0);
    for(int i = 0; i < data.size(); ++i){
        auto c0 = convol_c2r(vec_sc[i],w);
        for(int j = i; j < data.size(); ++j){
            auto c1 = convol_c2r(vec_sc[j],w);
            covar.push_back(inner_product(c0,c1,GridVol)*GridVol/datasize[j]/datasize[i] -1);
            delete[] c1;
        }
        delete[] c0;
    }
    std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
    for(auto x : covar) std::cout << x << ", "; std::cout << std::endl;
/*
    std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
    auto winpk = window_Pk(R,0);
    for(int i = 0; i < data.size(); ++i){
        for(int j = i; j < data.size(); ++j){
            auto pk_plus = densityCovarianceArray(vec_sc[i],vec_sc[j]);
            covar2.push_back(covar_CombinewithKernel(pk_plus,winpk));
            delete[] pk_plus;
        }
    }
    std::chrono::steady_clock::time_point time3 = std::chrono::steady_clock::now();
    for(auto x : covar2) std::cout << x << ", "; std::cout << std::endl;
    std::cout << "MRACS x-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << "[ms]" << std::endl;
    std::cout << "MRACS k-Statistics t = " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << "[ms]" << std::endl;
*/
    double* covariance = new double[16];
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            covariance[i * 4 + j] = covar[(i <= j) ? i * 4 - i*(i+1)/2 + j : j * 4 - j*(j+1)/2 + i];
/*    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            std::cout << covariance[i * 4 + j] << ", ";
        }
        std::cout << std::endl;
    }
*/
    double* eigen = new double[4];
    double* cov = new double[(4+1)*4/2];
    for(int i = 0; i < covar.size(); ++i) cov[i] = covar[i];

    auto info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',4,covariance,4,eigen);
    std::cout << "info: " << info << std::endl;
    for(int i = 0; i < 4; ++i) std::cout << eigen[i] << ", "; std::cout << std::endl;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j)
           std::cout << covariance[i * 4 + j] << ", ";
        std::cout << std::endl;
    }
}

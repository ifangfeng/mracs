// covariance of different halo bin 
#include"mracs.h"

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_uni = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(1,7);
    force_kernel_type(2);

    double GSR {2.1}; // Gaussian smoothing radius

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);
    
    //std::vector<double> vec_lambda_opt {4.5, 7.5, 9.25, 16.75, 17.75, 24.0, 23.75};
    double l_th {9.25};
    auto envi =  web_classify(cxx,hl,l_th);

    delete[] w_gs;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ----reconstruct and check-------
    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    const int Mbin{4};
    auto vpts_M = halo_mass_split(hl,Mbin);
    auto vpts_ME= halo_envi_mass_multi_split(envi, hl, Mbin);

    auto s_dm  = sfc(dm);
    auto s_hl  = sfc(hl);
    auto s_hl_u= sfc(hl_uni);

    auto sc_hl_uni = sfc_r2c(s_hl_u,false);
    auto sc_hl = sfc_r2c(s_hl,false);
    auto sc_dm = sfc_r2c(s_dm,false);

    auto sc_rc_M = optimal_reconstruct(sc_dm,vpts_M,Radius,true);
    auto s_rc_M  = sc_back(sc_rc_M,false);
    auto sc_rc_ME = optimal_reconstruct(sc_dm,vpts_ME,Radius,true);
    auto s_rc_ME  = sc_back(sc_rc_ME,false);

    std::vector<fftw_complex*> vec_sc {sc_dm,sc_hl_uni,sc_hl,sc_rc_M,sc_rc_ME};
    std::vector<double> norm;
    auto wpk_deltafunc = new double [(GridLen+1)*(GridLen+1)*(GridLen+1)];
    for(size_t i = 0; i < (GridLen+1)*(GridLen+1)*(GridLen+1); ++i) wpk_deltafunc[i] = 1;
    for(auto x : vec_sc){
        norm.push_back(covar_CombinewithKernel(densityVarianceArray(x),wpk_deltafunc,false));
    }

    for(int i = 0; i < GridLen/2 +1; ++i) std::cout << i* TWOPI / SimBoxL << ", ";
    std::cout << "Pk: dm, hl_uni, hl, hl_rc_M, hl_rc_ME\n";
    for(int i = 0; i < vec_sc.size(); ++i)
    {
        auto pk = densityPowerFFT(vec_sc[i]);
        for(int k = 0; k < GridLen/2 + 1; ++k){
            std::cout << pk[k] << ", ";
        }std::cout << std::endl;
        delete[] pk;
    }
    std::cout << "norm(aka. var): \n";
    for(auto x : norm) std::cout << x << ", "; std::cout << std::endl;

    force_kernel_type(0);
    auto vec_r = log_scale_generator(4,140,100,true);
    std::vector<double> xi_dm,xi_hl_u,xi_hl,xi_rc_M,xi_rc_ME;
    for(auto r : vec_r){
        auto w = wfc(r,0);
        auto c_dm = convol_c2r(sc_dm,w);
        auto c_hl_u = convol_c2r(sc_hl_uni,w);
        auto c_hl = convol_c2r(sc_hl,w);
        auto c_rc_M = convol_c2r(sc_rc_M,w);
        auto c_rc_ME= convol_c2r(sc_rc_ME,w);

        xi_dm.push_back(inner_product(s_dm,c_dm,GridVol)/pow(array_sum(s_dm,GridVol),2)*GridVol - 1);
        xi_hl_u.push_back(inner_product(s_hl_u,c_hl_u,GridVol)/pow(array_sum(s_hl_u,GridVol),2)*GridVol - 1);
        xi_hl.push_back(inner_product(s_hl,c_hl,GridVol)/pow(array_sum(s_hl,GridVol),2)*GridVol - 1);
        xi_rc_M.push_back(inner_product(s_rc_M,c_rc_M,GridVol)/pow(array_sum(s_rc_M,GridVol),2)*GridVol - 1);
        xi_rc_ME.push_back(inner_product(s_rc_ME,c_rc_ME,GridVol)/pow(array_sum(s_rc_ME,GridVol),2)*GridVol - 1);

        delete[] w;
        delete[] c_dm;
        delete[] c_hl_u;
        delete[] c_hl;
        delete[] c_rc_M;
        delete[] c_rc_ME;
    }
    std::cout << "correlation:\n";
    for(auto x : vec_r) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_dm) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_hl_u) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_hl) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_rc_M) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_rc_ME) std::cout << x << ", "; std::cout << std::endl;

}





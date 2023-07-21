// Hinfo  (double: x, y, z, Mass; int: M, C, E, S)
#include"mracs.h"

int main(int argc, char** argv){
    read_parameter();
    
    const int Ebin{4};
    int Mbin{1}, Cbin{1}, Sbin{1};
    
    if(argc==4){
        for(int i = 1; i < 4; ++i){
            int n = std::stoi(argv[i]);
            if(n<1 || n >128){
                std::cout << "input error, abort\n";
                return 0;
            }
        }
    }
    int Multi_size {Mbin * Cbin * Ebin * Sbin};
    if(Multi_size > 1024) {
        std::cout << "parameter space reaching the limit\n";
        std::terminate();
    }
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto h6 = read_in_Halo_6vector("/data0/MDPL2/halo_Mcut_slice.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    auto envi = envi_vector_readin(ifname,h6.size());

    // ----- classify halo catalogue then assign to Hinfo -------
    std::vector<Hinfo> hi(h6.size());
    #pragma omp parallel for
    for(size_t i = 0; i < h6.size(); ++i){
        hi[i] = {h6[i].x,h6[i].y,h6[i].z,h6[i].Mass,0,0,envi[i],0};
    }
    std::vector<double> vec_mass(h6.size());
    std::vector<double> vec_cont(h6.size());
    std::vector<double> vec_spin(h6.size());
    #pragma omp parallel for
    for(size_t i = 0; i < h6.size(); ++i){
        vec_mass[i] = h6[i].Mass;
        vec_cont[i] = h6[i].Concentration;
        vec_spin[i] = h6[i].Spin;
    }
    auto node_M = nodes_of_proto_sort(vec_mass, Mbin);
    auto node_C = nodes_of_proto_sort(vec_cont, Cbin);
    auto node_S = nodes_of_proto_sort(vec_spin, Sbin);

    for(size_t i = 0; i < h6.size(); ++i){
        hi[i].Mi = classify_index(node_M,h6[i].Mass);
        hi[i].Ci = classify_index(node_C,h6[i].Concentration);
        hi[i].Si = classify_index(node_S,h6[i].Spin);
    }

    std::vector<std::vector<Particle>*> vpts_mul;
    for(int i = 0; i < Multi_size; ++i) vpts_mul.push_back(new std::vector<Particle>);

    for(size_t i = 0; auto x : h6){
        vpts_mul[classify_index(node_M,x.Mass) * Cbin * Ebin * Sbin + classify_index(node_C,x.Concentration) * Ebin * Sbin + envi[i] * Sbin + 
        classify_index(node_S,x.Spin)]->push_back({x.x,x.y,x.z,x.Mass});
        ++i;
    }



    // ------vpts_uni-------
    std::vector<std::vector<Particle>*> vpts_uni;
    for(auto x : vpts) vpts_uni.push_back(new std::vector<Particle>);
    for(int i = 0; i < vpts.size(); ++i) for(auto x : *vpts[i]) vpts_uni[i]->push_back({x.x,x.y,x.z,1.});

    // ----reconstruct and check-------
    auto sc_hl_uni = sfc_r2c(sfc(hl_uni),true);
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_dm = sfc_r2c(sfc(dm),true);

    auto sc_rc_uni = optimal_reconstruct(dm,vpts_uni,Radius,true);
    auto sc_rc = optimal_reconstruct(dm,vpts,Radius,true);
    
    auto wpk = window_Pk(Radius,0);
    
    auto cc_uni = correlation_coefficients(sc_dm,sc_hl_uni,wpk);
    auto cc0 = correlation_coefficients(sc_dm,sc_rc_uni,wpk);
    auto cc1 = correlation_coefficients(sc_dm,sc_hl,wpk);
    auto cc2 = correlation_coefficients(sc_dm,sc_rc,wpk);
    
    std::cout << "Cross-correlation coefficient:\n";
    std::cout << "[number picture] default r_n: " << cc_uni << ", E= " << sqrt(1-pow(cc_uni,2)) << "\n";
    std::cout << "-----------RCST----------r_n: " << cc0    << ", E= " << sqrt(1-pow(cc0,2))    << "\n";
    std::cout << "[mass picture]   default r_m: " << cc1    << ", E= " << sqrt(1-pow(cc1,2))    << "\n";
    std::cout << "-----------RCST----------r_m: " << cc2    << ", E= " << sqrt(1-pow(cc2,2))    << "\n";
    


    
}





// covariance of vector of particle catalogue smoothed with radius R (kenel type specified in global variable "KernelFunc").
// if the vector in size n, then the returned vector has length (n+1)n/2, cov[i,j]==<vpts[i],vpts[j]>
std::vector<double> covar_of_data_vector_tmp(std::vector<std::vector<Particle>*> vpts, double R)
{
    std::vector<double> cov;
    //auto cov = new std::vector<double>;

    auto wpk = window_Pk(R,0);

    std::vector<fftw_complex*> vec_sc;

    for(auto x : vpts) 
        vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk,true));
    
    return cov;
}

// covariance of vector of particle catalogue. if the vector in size n,
// then the returned vector has length (n+1)n/2, cov[i,j]==<vpts[i],vpts[j]>
//std::vector<double> covar_of_data_vector_Multi_radii(std::vector<std::vector<Particle>*> vpts, std::vector<double>){
//
//}

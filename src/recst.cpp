// covariance of different halo bin 
#include"mracs.h"


int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};

    auto vpts = halo_envi_match_and_split(ifname,hl,dm);
    auto cov  = covar_of_data_vector(vpts,Radius);    
    auto solve = optimal_weight_solver(cov,4,true);
    std::cout << "solve:\n"; for(auto x : solve) std::cout << x << ", "; std::cout << "\n";

    
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
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk));
    
    return cov;
}

// covariance of vector of particle catalogue. if the vector in size n,
// then the returned vector has length (n+1)n/2, cov[i,j]==<vpts[i],vpts[j]>
//std::vector<double> covar_of_data_vector_Multi_radii(std::vector<std::vector<Particle>*> vpts, std::vector<double>){
//
//}

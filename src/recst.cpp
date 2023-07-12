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

    // --------reconstruct----------
    auto s = sfc(hl);
    std::vector<double> norm;
    double sum = static_cast<double>(hl.size());
    for(int i = 1; i < vpts.size(); ++i) norm.push_back(vpts[i]->size()/sum);
    for(size_t i = 0; i < hl.size(); ++i) hl[i].weight = 0;

    
}

// return the fourier of scaling function coefficients of optimal reconstructed halo catalogues, vpts is the
// vector of dm and splited halo cataloges, in envi-split case: {dm,vd,st,fl,kt}. R specify the smmothing scale
// which decide the reconstruct coeefficient of each halo component (weight vector), after solving weight
// vector we then reconstruct the halo fields with optimal weight and return as fourier of sfc coefficients
fftw_complex* optimal_reconstruct(std::vector<std::vector<Particle>*> vpts, double R)
{
    // ------covariance of each component------
    std::vector<double> cov;

    auto wpk = window_Pk(R,0);

    std::vector<fftw_complex*> vec_sc;

    for(auto x : vpts) 
        vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk));

    // ------solving optimal weight vector------
    int dim{vpts.size() - 1}; 
    auto solve =  optimal_weight_solver(cov,dim,true);

    // -----------------reconstruct--------------
    std::vector<double> weight, norm;
    for(int i = 1; i < solve.size(); ++i) weight.push_back(solve[i]);
    for(int i = 1; i < vpts.size(); ++i)  norm.push_back(static_cast<double>(vpts[i]->size()));

    double sum{0};
    for(auto x : norm) sum += x;
    for(int i = 0; i < weight.size(); ++i) weight[i] /= norm[i] / sum;

    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    #pragma omp parallel for
    for(size_t l = 0; l < GridLen * GridLen * (GridLen/2 + 1); ++l) { 
        sc[l][0] = 0; sc[l][1] = 0;}

    for(int i = 0; i < vec_sc.size() - 1; ++i)
        for(size_t l = 0; l < GridLen * GridLen * (GridLen/2 + 1); ++l) 
        {
            sc[l][0] += vec_sc[i+1][l][0];
            sc[l][1] += vec_sc[i+1][l][1];
        }

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

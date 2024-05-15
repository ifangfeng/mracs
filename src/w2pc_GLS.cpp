#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto pk = read_in_double("/home/feng/fac/data/Pk_Zeldovich_Planck13.bin");

    auto vec_r = linear_scale_generator(1,150,40,true);
    std::vector<double> vec_A{3,8,15}, xi_ms, xi_th;

    // --------- w2pc measure ----------
    auto s = sfc(p);
    auto sc = sfc_r2c(s,false);

    force_kernel_type(5);
    for(auto A : vec_A) {
        for(auto r : vec_r) {
            auto w = wfc(r,A);
            auto c = convol_c2r(sc,w);
            double dd = inner_product(c,s,GridVol);
            xi_ms.push_back(dd * GridVol/pow(p.size(), 2) - 1);
            delete[] w;
            delete[] c;
        }
    }
    
    for(auto r : vec_r) std::cout << r << ", "; std::cout  << std::endl;

    for(int i = 0; i < vec_A.size(); ++i) {
        for(int j = 0; j < vec_r.size(); ++j){
            std::cout << xi_ms[i*vec_r.size() + j] << ", "; 
        }
        std::cout << std::endl;
    }

    // --------- w2pc theory ----------
    auto vec_r2 = linear_scale_generator(1,150,150,true);
    for(auto A : vec_A){
        for(auto r : vec_r2){
            auto win_GLS = win_theory("GLS",r,A);
            double sum{0};
            #pragma omp parallel for reduction (+:sum)
            for(int i = 0; i < pk.size(); ++i){
                sum += pk[i]*win_GLS[i]*pow((i+1)*1e-4,2);
            }
            xi_th.push_back((sum*1e-4)*4*M_PI/pow(TWOPI,3));
        } 
    }

    for(auto A : vec_A) 
        std::cout << "Scale A:" << A << ", "; std::cout  << std::endl;

    for(auto r : vec_r2) std::cout << r << ", "; std::cout  << std::endl;

    for(int i = 0; i < vec_A.size(); ++i){
        for(int j = 0; j < vec_r2.size(); ++j){
            std::cout << xi_th[i*vec_r2.size()+j] << ", ";
        }std::cout << std::endl;
    }

}
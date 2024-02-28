#include"mracs.h"


int main()
{
    read_parameter();
    auto pk = read_in_double("/home/feng/fac/data/Pk_Planck13.bin");

    auto vec_r = linear_scale_generator(1,150,40,true);
    std::vector<double> vec_R {3,8,15}, var;
    for(auto R : vec_R){
        auto win_R = win_theory("GDW",R);
        std::vector<double> w_tmp(pk.size());
        #pragma omp parallel for
        for(int i = 0; i < pk.size(); ++i)
            w_tmp[i] = pk[i] * win_R[i] * win_R[i];
        for(auto r : vec_r){
            auto win_r = win_theory("sphere",r);
            double sum{0};
            #pragma omp parallel for reduction (+:sum)
            for(int i = 0; i < pk.size(); ++i){
                sum += w_tmp[i]*win_r[i]*win_r[i]*pow((i+1)*1e-4,2);
            }
            var.push_back((sum*1e-4)*4*M_PI/pow(TWOPI,3));
        } 
    }

    for(auto R : vec_R) 
        std::cout << "smooth R:" << R << ", "; std::cout  << std::endl;
    for(int i = 0; i < vec_R.size(); ++i){
        for(int j = 0; j < vec_r.size(); ++j){
            std::cout << var[i*vec_r.size()+j] << ", ";
        }std::cout << std::endl;
    }

}



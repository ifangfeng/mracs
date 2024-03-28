#include"mracs.h"


int main()
{
    read_parameter();
    auto pk = read_in_double("/home/feng/fac/data/Pk_Zeldovich_Planck13.bin");

    auto vec_R = linear_scale_generator(1,200,100,true);
    std::vector<double> var;
    for(auto R : vec_R){
        auto win_R = win_theory("GDW",R,0);
        double sum{0};
        #pragma omp parallel for reduction (+:sum)
        for(int i = 0; i < pk.size(); ++i){
            sum += pk[i]*win_R[i]*win_R[i]*pow((i+1)*1e-4,2);
        }
        var.push_back((sum*1e-4)*4*M_PI/pow(TWOPI,3));
    }

    for(auto x : var) std::cout << x << ", "; std::cout << std::endl;

}



#include"mracs.h"

std::vector<double> win_theory(std::string kernel, double R);

int main()
{
    read_parameter();
    auto p = read_in_TNG_3vector("/data0/BigMDPL/dm_particles_snap_079_position.bin");
    auto pk = read_in_double("/home/feng/fac/data/Pk_Planck15.bin");
    auto win = win_theory("shell",10)
    
    for(int i = 0; i < )

}

// sampling the window fourier from [1e-4,10] with delta_k = 1e-4
// combine with Pk_theory value to calculate the theoretical second order statistics
// kernel can be "shell" "Gaussian","sphere","GDW"
std::vector<double> win_theory(std::string kernel, double R){
    const double k_min = 1e-4; //radian frequency
    const double k_max = 10;
    const double dk = 1e-4;
    const int klen = 1e5;

    std::vector<double> win(klen);
    if(kernel == "shell") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = sin(phase)/(phase);
        }
    }
    else if(kernel == "sphere") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = 3*(sin(phase)-phase*cos(phase))/(pow(phase,3));
        }
    }
    else if(kernel == "Gaussian") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = pow(1/M_E,phase*phase/2);
        }
    }
    else if(kernel == "GDW") {
        const double norm = pow(2,7./4)/sqrt(15)*pow(TWOPI,3./4)*pow(R,3./2);
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = sin(phase)/(phase);
        }
    }

    return win;
}
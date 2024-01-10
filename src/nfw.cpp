#include"mracs.h"

double M2R(double M){
    const double rhoHalo {2.83386e13}; //M_star*h^-1/[Mpc*h^-1]^3
    return pow(M/(4./3*M_PI*rhoHalo),1./3);
}

void zero_check(double x){
    if(!x) std::cout << "zero\n";
    else std::cout << "non-zero\n"; 
}

int main(){
    read_parameter();


    std::vector<double> haloMass {1e13, 1e14, 1e15};
    std::vector<double> haloConce {6.9, 5.7, 4.7};

    const double deltaXi {1./GridLen};
    for(int i = 0; i < haloMass.size(); ++i){

        double r_h = M2R(haloMass[i]);
        double r_s = r_h/haloConce[i];

        std::cout << "r_h: " << r_h << ", r_s:" << r_s << std::endl;

        for(int j = 0; j < GridLen; ++j){
            std::cout << NFW_window_norm(r_h * GridLen/SimBoxL, r_s * GridLen/SimBoxL,10, 0,0,j*deltaXi) << ", ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "k_Nyq=" << M_PI/(SimBoxL/GridLen) << std::endl;
}


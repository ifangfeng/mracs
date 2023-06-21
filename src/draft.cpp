#include"mracs.h"

int main(){
    read_parameter();

    double array[9]; for(int i = 0; i < 9; ++i) array[i] = i+1; array[8] = 10;
    double array_b[6]{0};array_b[0] = 1;array_b[3]=1;
    int ipiv[3];

    std::cout << "ipiv0:\n";for(int i = 0; i < 3; ++i) std::cout << ipiv[i] << ", ";std::cout << std::endl;

    //auto info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, 3,2,array,3,ipiv,array_b,2);
    auto info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,3,3,array,3,ipiv);
    std::cout << "ipiv1:\n";for(int i = 0; i < 3; ++i) std::cout << ipiv[i] << ", ";std::cout << std::endl;
    auto info2= LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',3,2,array,3,ipiv,array_b,2);

    std::cout << "ipiv2:\n";for(int i = 0; i < 3; ++i) std::cout << ipiv[i] << ", ";std::cout << std::endl;
    std::cout << "info: " << info << ", " << info2 << std::endl;
    std::cout << "array_b: "; for(int i = 0; i < 6; ++i) std::cout << array_b[i] << ", "; std::cout << std::endl;

    
}
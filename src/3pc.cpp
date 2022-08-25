#include"mracs.h"

using namespace std;

double WindowFunction_Dual_Ring1(double R, double theta, double ki, double kj, double kk);

int main()
{
    int a = 3;

    int L = 128;
    size_t N = 128;
    auto c1 = a & (L-1);
    auto d1 = a & (N-1);

    std::cout << c1 << ", " << d1 << ", " << std::endl;
    for(int i = 0; i < 100; ++i)
    {
        std::cout << i/100. << ", ";
    }
    std::cout << "\n";
}


//dual-ring window function: J_0(2Pi*kRsin(a)sin(b))*cos(2Pi*kRcos(a)cos(b))
double WindowFunction_Dual_Ring1(double R, double theta, double ki, double kj, double kk)
{
    return (std::cyl_bessel_j(0,TWOPI*sin(theta)*R*sqrt(ki*ki+kj*kj)))*cos(TWOPI * R * cos(theta) * kk);
}
#include"mracs.h"

int main()
{
    int a = 3;

    int L = 128;
    size_t N = 128;
    auto c1 = a & (L-1);
    auto d1 = a & (N-1);

    std::cout << c1 << ", " << d1 << ", " << std::endl;
}
#include"mracs.h"

std::vector<double> kd_tree(std::vector<Particle> p)
{
    const int len = log2(p.size());
    const int len_2d = (1 + len) / 2;
    std::vector<double> a;
    double temp{0};

    for(int i = 0; i < len; ++i)
        for(int j = 0; j < 1 << i; ++j)
        {

            for(size_t n = j * (1 << i); n < p.size(); ++n) temp += p[n].x;
            a.push_back(temp/2);
        }
}


std::vector<int64_t> count_with_tree(std::vector<Particle> p, std::vector<Point> pt)
{
    std::vector<int64_t> result;

    for(size_t i = 0; i < pt.size(); ++i)
    {

    }
}
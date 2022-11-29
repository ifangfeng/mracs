#include<iostream>
#include<string>
#include<vector>
#include<cmath>

struct tree
{
    unsigned l;
};

struct Point
{
    double x;
    double y;
    double z;
};


int main()
{
    std::vector<Point> pt;
    unsigned N = log2(pt.size());
    uint64_t ID;
    
    for(int i = 0; i < sizeof(ID)*8; ++i)
        if(ID & (1 << i))

}
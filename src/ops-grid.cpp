// ********************************
// One Point Statistics--FFTbase
// ********************************

#include"mracs.h"
#include"kdtree.hpp"
using namespace std;

int main()
{
    read_parameter();
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    
    force_base_type(0,1);
    force_kernel_type(1);
    auto s = sfc(p1);
    auto sc = sfc_r2c(s);
    auto w = wft(Radius,0);
    auto c = convol_c2r(sc,w);      // FFT-base density filed on grids

    auto w1 = wfc(Radius,0);
    auto cs = convol_c2r(sc,w1);
    auto c1 = prj_grid(cs);         // DWT-base density filed on grids

    
}
// *******************************************************
// print halo environment information at its nearst grid point
// *******************************************************
#include"mracs.h"
using namespace std;

int main(){
    read_parameter();
    int GSR {2}; // Gaussian smoothing radius
    std::string ofname {"output/envi_J" + std::to_string(Resolution)+ "_GSR"+std::to_string(GSR) +"_halo_Mcut2e12.txt"};
    std::ofstream envirol {ofname};

    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");


    force_base_type(0,1);
    force_kernel_type(2);
    auto sc = sfc_r2c(sfc(p1),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p2);

    for(auto x : env)
        envirol << x << ", ";
    std::cout << std::endl;

    int64_t voids{0}, sheets{0}, filaments{0}, knots{0};
    for(auto x : env){
        if (x == 0) voids++;
        else if (x == 1) sheets++;
        else if (x == 2) filaments++;
        else knots++;
    }

    double sum = voids + sheets + filaments + knots;
    std::cout << "number of halo : " << p2.size() << std::endl;
    std::cout << "sum: " << sum << std::endl;
    std::cout << "voids: " << voids/sum << std::endl;
    std::cout << "sheets: " << sheets/sum << std::endl;
    std::cout << "filaments: " << filaments/sum << std::endl;
    std::cout << "knots: " << knots/sum << std::endl;
    

}
// *******************************************************
// print halo environment information at its nearst grid point
// *******************************************************
#include"mracs.h"
using namespace std;

int main(){
    read_parameter();
    int GSR {5}; // Gaussian smoothing radius
    std::string ofname {"output/envi_J" + std::to_string(Resolution)+ "_GSR"+std::to_string(GSR) +"_halo_Mcut2e12_grid.txt"};
    std::ofstream envirol {ofname};

    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    //auto p0 = generate_random_particle(1024,SimBoxL,0);

    force_base_type(0,1);
    force_kernel_type(2);
    auto s = sfc(p1);
    auto sc = sfc_r2c(s);
    auto w = wft(GSR, 0);

    auto begin = std::chrono::steady_clock::now();
    auto cxx = tidal_tensor(sc, w);
    auto end1 = std::chrono::steady_clock::now();
    auto env = web_classify_to_grid(cxx);
    auto end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference tidal    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin).count() << "[ms]" << std::endl;
    std::cout << "Time difference classify    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << "[ms]" << std::endl;

    auto begin3 = std::chrono::steady_clock::now();
    for(auto x : env)
        envirol << x << ", ";
    auto end3 = std::chrono::steady_clock::now();
    std::cout << "Time difference write    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count() << "[ms]" << std::endl;


    int64_t voids{0}, sheets{0}, filaments{0}, knots{0};
    for(auto x : env){
        if (x == 0) voids++;
        else if (x == 1) sheets++;
        else if (x == 2) filaments++;
        else knots++;
    }

    double sum = voids + sheets + filaments + knots;
    std::cout << "GridVol: " << GridVol << std::endl;
    std::cout << "sum: " << sum << std::endl;
    std::cout << "voids: " << voids/sum << std::endl;
    std::cout << "sheets: " << sheets/sum << std::endl;
    std::cout << "filaments: " << filaments/sum << std::endl;
    std::cout << "knots: " << knots/sum << std::endl;
    

}
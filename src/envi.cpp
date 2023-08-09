// *******************************************************
// print halo environment information at its nearst grid point
// *******************************************************
#include"mracs.h"
using namespace std;

int main(){
    read_parameter();
    double GSR {2.1}; // Gaussian smoothing radius
   
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    force_kernel_type(2);
    auto sc = sfc_r2c(sfc(p1),true);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p2,0);

    std::vector<double> vd,st,fl,kt;
    std::vector<double> vd2,st2,fl2,kt2;
    for(int i = 0; i < 20; ++i){
        auto envGrid = web_classify_to_grid(cxx,10*i/20.);
        auto env = web_classify(cxx,p2,10*i/20.);
        int64_t voids{0}, sheets{0}, filaments{0}, knots{0};
        int64_t voids2{0}, sheets2{0}, filaments2{0}, knots2{0};
        for(auto x : envGrid){
            if (x == 0) voids++;
            else if (x == 1) sheets++;
            else if (x == 2) filaments++;
            else knots++;
        }
        for(auto x : env){
            if (x == 0) voids2++;
            else if (x == 1) sheets2++;
            else if (x == 2) filaments2++;
            else knots2++;
        }
        std::vector<int>().swap(envGrid);
        std::vector<int>().swap(env);
        double sum = voids + sheets + filaments + knots;
        double sum2 = voids2 + sheets2 + filaments2 + knots2;
        vd.push_back(voids/sum);
        st.push_back(sheets/sum);
        fl.push_back(filaments/sum);
        kt.push_back(knots/sum);
        vd2.push_back(voids2/sum2);
        st2.push_back(sheets2/sum2);
        fl2.push_back(filaments2/sum2);
        kt2.push_back(knots2/sum2);
        std::cout << i << "\n";
    }
    std::ofstream ofn_vff_vd{"output/VFF_vd.txt"},ofn_HF_vd{"output/HF_vd.txt"};
    std::ofstream ofn_vff_st{"output/VFF_st.txt"},ofn_HF_st{"output/HF_st.txt"};
    std::ofstream ofn_vff_fl{"output/VFF_fl.txt"},ofn_HF_fl{"output/HF_fl.txt"};
    std::ofstream ofn_vff_kt{"output/VFF_kt.txt"},ofn_HF_kt{"output/HF_kt.txt"};

    for(auto x : vd) ofn_vff_vd << x << " ";
    for(auto x : st) ofn_vff_st << x << " ";
    for(auto x : fl) ofn_vff_fl << x << " ";
    for(auto x : kt) ofn_vff_kt << x << " ";
    for(auto x : vd2) ofn_HF_vd << x << " ";
    for(auto x : st2) ofn_HF_st << x << " ";
    for(auto x : fl2) ofn_HF_fl << x << " ";
    for(auto x : kt2) ofn_HF_kt << x << " ";

}
#include"mracs.h"

void fraction_growing(std::vector<double>& vec_frac, std::vector<int>& env){

    double norm = env.size();
    size_t vd{0},st{0},fl{0},kt{0};
    for(auto x : env){
        if(x == 0) vd++;
        else if(x == 1) st++;
        else if(x == 2) fl++;
        else if(x == 3) kt++;
        else {std::cout << "Func[fraction_growing]: env input error! abort\n"; std::terminate();}
    }
    
    vec_frac.push_back(vd / norm);
    vec_frac.push_back(st / norm);
    vec_frac.push_back(fl / norm);
    vec_frac.push_back(kt / norm);
}

int main(){
    read_parameter();

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub05.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::ofstream ofs_n {"output/rcnf.txt"};
    std::ofstream ofs_v {"output/rcvf.txt"};

    double GSR {1}; // Gaussian smoothing radius

    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);

    delete[] w_gs;
    fftw_free(sc);

    std::vector<double> VF, NF;
    auto vec_lth = log_scale_generator(0.1,20,40,true);

    for(auto l_th : vec_lth) {
        auto env =  web_classify(cxx,hl,l_th);
        auto env_grid = web_classify_to_grid(cxx,l_th);
        fraction_growing(NF,env);
        fraction_growing(VF,env_grid);
        std::vector<int>().swap(env);
        std::vector<int>().swap(env_grid);
    }

    for(auto x : VF) ofs_v << x << " "; ofs_v << std::endl;
    for(auto x : NF) ofs_n << x << " "; ofs_n << std::endl;
}


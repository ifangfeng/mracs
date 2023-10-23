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
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    double GSR {1}; // Gaussian smoothing radius
    double THR {30};

    std::ofstream ofs_omn {"output/rcom_THR30_n.txt"};
    std::ofstream ofs_omm {"output/rcom_THR30_m.txt"};
    std::ofstream ofs_NF  {"output/rcom_THR30_NF.txt"};
    std::ofstream ofs_Msum{"output/rcom_THR30_binMass.txt"};
    std::ofstream ofs_ME  {"output/rcom_THR30_dataME.txt"};

    const int Mbin{4};
    double lth_opt {17.5};
    // ------environment sticker---------
    force_resoluton_J(10);
    force_base_type(0,1);
    force_kernel_type(2);

    auto sc = sfc_r2c(sfc(dm),true);
    auto w_gs = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w_gs);

    auto env_ME = web_classify(cxx,hl,lth_opt);
    auto vec_lth = log_scale_generator(0.1,20,40,true);
    std::vector<std::vector<int>> vec_env;
    for(auto l_th : vec_lth) 
        vec_env.push_back(web_classify(cxx,hl,l_th));
    
    delete[] w_gs;
    fftw_free(sc);
    for(int i = 0; i < 6; ++i) delete[] cxx[i];delete[] cxx;

    // ------cross-correlation of different Lambda_th---------
    force_resoluton_J(9);
    force_base_type(1,4);
    force_kernel_type(1);

    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto wpk = window_Pk(THR,0);

    // -----line one---
    std::vector<double> W_M, NF, Msum, W_N, data_ME;
    for(auto env : vec_env) {
        auto vpts_n = halo_envi_match_and_split(env,hl_n);
        auto sol_n = optimal_solution_lean(sc_dm,vpts_n,wpk,false);

        auto vpts_m = halo_envi_match_and_split(env,hl);
        auto sol_m = optimal_solution_lean(sc_dm,vpts_m,wpk,false);
        for(int i = 1; i < sol_n.size(); ++i) W_N.push_back(sol_n[i]);
        for(int i = 1; i < sol_m.size(); ++i) W_M.push_back(sol_m[i]);

        for(auto x : vpts_m) {
            double sum{0}; for(auto pt : *x) sum += pt.weight;
            Msum.push_back(sum);
        }
        fraction_growing(NF,env);
        std::vector<int>().swap(env);
    }
    // ------Mbin=4, Ebin=4------
    auto vpts = halo_envi_mass_multi_split(env_ME,hl,Mbin);
    auto sol = optimal_solution_lean(sc_dm,vpts,wpk,false);
    for(int i = 1; i < sol.size(); ++i) data_ME.push_back(sol[i]);
    for(auto x : vpts) data_ME.push_back(x->size()/ static_cast<double>(hl.size()));
    for(auto x : vpts) {
        double sum{0}; for(auto pt : *x) sum += pt.weight;
        data_ME.push_back(sum);
    }
    for(int i = 0; auto x : data_ME){
        std::cout << x << ", ";
        if((i+1)%4==0) std::cout << std::endl;
        ++i;
    }

    // ------output------
    for(auto x : W_M)    ofs_omm  << x << " "; ofs_omm  << std::endl;
    for(auto x : W_N)    ofs_omn  << x << " "; ofs_omn  << std::endl;
    for(auto x : Msum)    ofs_Msum << x << " "; ofs_Msum << std::endl;
    for(auto x : NF)      ofs_NF   << x << " "; ofs_NF   << std::endl;
    for(auto x : data_ME) ofs_ME   << x << " "; ofs_ME   << std::endl;

}
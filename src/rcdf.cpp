// density map strips
#include"mracs.h"


int main(int argc, char** argv){
    read_parameter();
    
    int Mbin{0};
    const int Ebin{4};

    if(argc==1){
        Mbin = 4;
    }
    else {
        int n = std::stoi(argv[1]);
        if ((n >= 1) && (n <= 128))
            Mbin = n;
        else {
            std::cout << "input error, abort\n";
            return 0;
        }
    }

    auto dm   = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl_n = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_m = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname      {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    std::string ofname_dm   {"output/rcdf_dm_THR"+RADII+".txt"};
    std::string ofname_hl_n {"output/rcdf_hl_n_THR"+RADII+".txt"};
    std::string ofname_hl_m {"output/rcdf_hl_m_THR"+RADII+".txt"};
    std::string ofname_rc_M {"output/rcdf_M_split_THR" + RADII + ".txt"};
    std::string ofname_rc_ME{"output/rcdf_ME_split_THR" + RADII + ".txt"};

    std::ofstream ofsdm{ofname_dm}, ofshl_n{ofname_hl_n}, ofshl_m{ofname_hl_m},
    ofsrc_M{ofname_rc_M}, ofsrc_ME{ofname_rc_ME};

    std::vector<std::vector<Particle>*> vpts_M, vpts_ME;

    vpts_M = halo_mass_split(hl_m, Mbin);
    vpts_ME = halo_envi_mass_multi_split(ifname, hl_m, Mbin);

    
    // ----reconstruct-------
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_hl_m = sfc_r2c(sfc(hl_m),true);
    auto sc_hl_n = sfc_r2c(sfc(hl_n),true);

    auto sc_rc_M = optimal_reconstruct(dm,vpts_M,Radius,true);
    auto sc_rc_ME = optimal_reconstruct(dm,vpts_ME,Radius,true);

    auto w = wfc(Radius,0);
    auto s_dm = convol_c2r(sc_dm, w);
    auto s_hl_m = convol_c2r(sc_hl_m, w);
    auto s_hl_n = convol_c2r(sc_hl_n, w);

    auto s_rc_M = convol_c2r(sc_rc_M, w);
    auto s_rc_ME = convol_c2r(sc_rc_ME, w);
    
    // -----sampling------
    double sum_dm{0},sum_hl_m{0},sum_hl_n{0},sum_rc_M{0},sum_rc_ME{0};
    #pragma omp parallel for reduction (+:sum_dm)
    for(size_t i = 0; i < GridVol; ++i) sum_dm += s_dm[i];
    #pragma omp parallel for reduction (+:sum_hl_m)
    for(size_t i = 0; i < GridVol; ++i) sum_hl_m += s_hl_m[i];
    #pragma omp parallel for reduction (+:sum_hl_n)
    for(size_t i = 0; i < GridVol; ++i) sum_hl_n += s_hl_n[i];
    
    #pragma omp parallel for reduction (+:sum_rc_M)
    for(size_t i = 0; i < GridVol; ++i) sum_rc_M += s_rc_M[i];
    #pragma omp parallel for reduction (+:sum_rc_ME)
    for(size_t i = 0; i < GridVol; ++i) sum_rc_ME += s_rc_ME[i];

    // ---slice---
    std::vector<Particle> p0;
    const int Nslice {21};
    int Ly{200},Lz{1000};
    int    Ny{Ly},Nz{Lz};
    auto y = linear_scale_generator(0,Ly,Ny,false);
    auto z = linear_scale_generator(0,Lz,Nz,false);
    for(int i = 0; i < Nslice; ++i){
        std::string nid = std::to_string(i);
        std::string ofname_dm   {"output/rcdf_dm_THR"+RADII+"_x"+nid+".txt"};
        std::string ofname_hl_n {"output/rcdf_hl_n_THR"+RADII+"_x"+nid+".txt"};
        std::string ofname_hl_m {"output/rcdf_hl_m_THR"+RADII+"_x"+nid+".txt"};
        std::string ofname_rc_M {"output/rcdf_M_split_THR" + RADII+"_x"+nid + ".txt"};
        std::string ofname_rc_ME{"output/rcdf_ME_split_THR" + RADII+"_x"+nid + ".txt"};

        std::ofstream ofsdm{ofname_dm}, ofshl_n{ofname_hl_n}, ofshl_m{ofname_hl_m},ofsrc_M{ofname_rc_M}, ofsrc_ME{ofname_rc_ME};
        for(auto j : y)
            for(auto k : z){
                p0.push_back({i*SimBoxL/Nslice,j,k,1.});
            }
    

    auto n_dm = project_value(s_dm,p0,false);
    auto n_hl_m = project_value(s_hl_m,p0,false);
    auto n_hl_n = project_value(s_hl_n,p0,false);
    auto n_rc_M = project_value(s_rc_M,p0,false);
    auto n_rc_ME = project_value(s_rc_ME,p0,false);

    std::vector<Particle>().swap(p0);
    for(size_t i = 0; i < Ny*Nz; ++i) ofsdm << n_dm[i] / sum_dm * GridVol -1 << " ";ofsdm.close();
    for(size_t i = 0; i < Ny*Nz; ++i) ofshl_m << n_hl_m[i] / sum_hl_m * GridVol -1 << " ";ofshl_m.close();
    for(size_t i = 0; i < Ny*Nz; ++i) ofshl_n << n_hl_n[i] / sum_hl_n * GridVol -1 << " ";ofshl_n.close();
    for(size_t i = 0; i < Ny*Nz; ++i) ofsrc_M << n_rc_M[i] / sum_rc_M * GridVol -1 << " ";ofsrc_M.close();
    for(size_t i = 0; i < Ny*Nz; ++i) ofsrc_ME << n_rc_ME[i] / sum_rc_ME * GridVol -1 << " ";ofsrc_ME.close();

    delete[] n_dm;
    delete[] n_hl_m;
    delete[] n_hl_n;
    delete[] n_rc_M;
    delete[] n_rc_ME;
    }
    
}






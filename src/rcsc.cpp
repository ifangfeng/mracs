// scatter of reconstructed halo vs. dark matter catalogue
#include"mracs.h"


int main(int argc, char** argv){
    read_parameter();
    
    
    if(argc>1){
        int n = std::stoi(argv[1]);
        if(n<1 || n >128){
            std::cout << "input error, abort\n";
            return 0;
        }
    }

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    std::string ofname_dm {"output/rcss_dm_THR"+RADII+".txt"};
    std::string ofname_hl {"output/rcss_hl_num_THR"+RADII+".txt"};
    std::string ofname_rc;

    std::vector<std::vector<Particle>*> vpts;

    if(argc == 1){
        ofname_rc = "output/rcss_envi_split_THR" + RADII + ".txt";
        vpts = halo_envi_match_and_split(ifname,hl);
    }
    else{
        int n = std::stoi(argv[1]);
        ofname_rc = "output/rcss_mass_split_THR" + RADII +"_N"+ std::to_string(n) + ".txt";
        vpts = halo_mass_split(hl,n);
    }
    std::ofstream ofsdm {ofname_dm}, ofshl {ofname_hl}, ofsrc {ofname_rc};

    // ----reconstruct-------
    hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto sc_hl = sfc_r2c(sfc(hl),true);
    auto sc_dm = sfc_r2c(sfc(dm),true);
    auto sc_rc = optimal_reconstruct(dm,vpts,Radius,true);

    auto w = wfc(Radius,0);
    auto s_hl = convol_c2r(sc_hl, w);
    auto s_dm = convol_c2r(sc_dm, w);
    auto s_rc = convol_c2r(sc_rc, w);
    
    // -----sampling------
    double sum_dm{0},sum_hl{0},sum_rc{0};
    #pragma omp parallel for reduction (+:sum_hl)
    for(size_t i = 0; i < GridVol; ++i) sum_hl += s_hl[i];
    #pragma omp parallel for reduction (+:sum_dm)
    for(size_t i = 0; i < GridVol; ++i) sum_dm += s_dm[i];
    #pragma omp parallel for reduction (+:sum_rc)
    for(size_t i = 0; i < GridVol; ++i) sum_rc += s_rc[i];

    size_t N{80};
    auto p0 = generate_random_particle(N,SimBoxL,0);

    auto n_dm = project_value(s_dm,p0,false);
    auto n_hl = project_value(s_hl,p0,false);
    auto n_rc = project_value(s_rc,p0,false);

    
    for(size_t i = 0; i < N*N*N; ++i) ofsdm << n_dm[i] / sum_dm * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofshl << n_hl[i] / sum_hl * GridVol << " ";
    for(size_t i = 0; i < N*N*N; ++i) ofsrc << n_rc[i] / sum_rc * GridVol << " ";
    
}






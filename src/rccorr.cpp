// Hinfo  (double: x, y, z, Mass; int: M, E, C, D)
#include"mracs.h"
//---Number of Particles: 283116474 of sub005
int main(int argc, char** argv){
    // ----initial-----
    const int Ebin{4};
    int Mbin{1}, Cbin{1}, Dbin{1};
    
    if(argc!=4){
        std::cout << "Three parameters are needed\n Each bin parameter should be lager than or equal to 1, \n MassBin ConcentratrionBin LocDensBin (M C D)\n";
        std::cout << " e.g.\n";
        std::cout << argv[0] << " 4 4 4\n"; return 0;
    }
    else {
        for(int i = 1; i < 4; ++i){
            int n = std::stoi(argv[i]);
            if(n<1 || n >128){
                std::cout << "input error, abort\n";
                return 0;
            }
        }
        Mbin = std::stoi(argv[1]);
        Cbin = std::stoi(argv[2]);
        Dbin = std::stoi(argv[3]);
    }
    int Multi_size {Mbin * Ebin * Cbin * Dbin};
    if(Multi_size > 1024) {
        std::cout << "parameter space reaching the limit\n";
        std::terminate();
    }
    read_parameter();
    // ------reading----
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto h6 = read_in_Halo_6vector("/data0/MDPL2/halo_Mcut_slice.bin");
    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    auto envi = envi_vector_readin(ifname,h6.size());

    std::vector<Particle> hl(h6.size()), hl_uni(h6.size());
    for(size_t i = 0; i < h6.size(); ++i) {
        hl[i] = {h6[i].x,h6[i].y,h6[i].z,h6[i].Mass};
        hl_uni[i] = {h6[i].x,h6[i].y,h6[i].z,1.};
    }

    // ----substitute locDensity for spin----
    auto s = sfc(hl);
    auto w = wfc(Radius,0);
    auto c = convol3d(s,w,true);
    auto nprj = project_value(c,hl,true);
    delete[] w;
    for(size_t i = 0; i < h6.size(); ++i){
        h6[i].Spin = nprj[i];
    }


    // ======================= classify and dispatch halo particles ================
    std::vector<double> vec_mass(h6.size()), vec_cntr(h6.size()), vec_spin(h6.size());
    #pragma omp parallel for
    for(size_t i = 0; i < h6.size(); ++i){
        vec_mass[i] = h6[i].Mass;
        vec_cntr[i] = h6[i].Concentration;
        vec_spin[i] = h6[i].Spin;
    }
    auto node_M = nodes_of_proto_sort(vec_mass, Mbin);
    auto node_C = nodes_of_proto_sort(vec_cntr, Cbin);
    auto node_S = nodes_of_proto_sort(vec_spin, Dbin);

    std::cout << "spin:\n";
    std::vector<std::vector<double>*> vv_spin;
    for(int i = 0; i < Dbin; ++i) vv_spin.push_back(new std::vector<double>);
    for(auto x : vec_spin) vv_spin[classify_index(node_S,x)]->push_back(x);
    for(auto x : vv_spin)
        print_min_max_and_size_double(*x);


    
    
    std::vector<std::vector<Particle>*> vpts_mass, vpts_envi, vpts_cntr, vpts_spin;
    std::vector<std::vector<Particle>*> vpts_coca, vpts_mult;
    for(int i = 0; i < Dbin; ++i) vpts_spin.push_back(new std::vector<Particle>);
    for(int i = 0; i < Ebin; ++i) vpts_envi.push_back(new std::vector<Particle>);
    for(int i = 0; i < Cbin; ++i) vpts_cntr.push_back(new std::vector<Particle>);
    for(int i = 0; i < Mbin; ++i) vpts_mass.push_back(new std::vector<Particle>);
    for(int i = 0; i < Mbin + Cbin + Ebin + Dbin; ++i) vpts_coca.push_back(new std::vector<Particle>);
    for(int i = 0; i < Mbin * Cbin * Ebin * Dbin; ++i) vpts_mult.push_back(new std::vector<Particle>);
    

    // ------------------- vpts_ M C E S------------------
    for(auto x : h6) vpts_mass[classify_index(node_M,x.Mass)]->push_back({x.x,x.y,x.z,x.Mass});
    for(auto x : h6) vpts_cntr[classify_index(node_C,x.Concentration)]->push_back({x.x,x.y,x.z,x.Mass});
    for(size_t i = 0; i < h6.size(); ++i) vpts_envi[envi[i]]->push_back({h6[i].x,h6[i].y,h6[i].z,h6[i].Mass});
    for(auto x : h6) vpts_spin[classify_index(node_S,x.Spin)]->push_back({x.x,x.y,x.z,x.Mass});

    // ------------------- multi-dimension ------------------
    for(size_t i = 0; auto x : h6){
        vpts_mult[classify_index(node_M,x.Mass) * Cbin * Ebin * Dbin + classify_index(node_C,x.Concentration) * Ebin * Dbin + envi[i] * Dbin + 
        classify_index(node_S,x.Spin)]->push_back({x.x,x.y,x.z,x.Mass});
        ++i;
    }
    // ------------------- concatenate ------------------
    for(auto x : h6) vpts_coca[classify_index(node_M,x.Mass)]->push_back({x.x,x.y,x.z,x.Mass});
    for(auto x : h6) vpts_coca[Mbin + classify_index(node_C,x.Concentration)]->push_back({x.x,x.y,x.z,x.Mass});
    for(size_t i = 0; i < h6.size(); ++i) vpts_coca[Mbin + Cbin + envi[i]]->push_back({h6[i].x,h6[i].y,h6[i].z,h6[i].Mass});
    for(auto x : h6) vpts_coca[Mbin + Cbin + Ebin + classify_index(node_S,x.Spin)]->push_back({x.x,x.y,x.z,x.Mass});

    std::cout << "spin mass:\n";
    for(auto x : vpts_spin) print_min_max_and_size(*x);

    // ----reconstruct and check-------
    auto s_hl_uni = sfc(hl_uni);
    auto s_hl = sfc(hl);
    auto s_dm = sfc(dm);

    auto sc_hl_uni = sfc_r2c(s_hl_uni,false);
    auto sc_hl = sfc_r2c(s_hl,false);
    auto sc_dm = sfc_r2c(s_dm,false);

    auto sc_rc_mass = optimal_reconstruct(dm,vpts_mass,Radius,true);
    auto sc_rc_cntr = optimal_reconstruct(dm,vpts_cntr,Radius,true);
    auto sc_rc_envi = optimal_reconstruct(dm,vpts_envi,Radius,true);
    auto sc_rc_spin = optimal_reconstruct(dm,vpts_spin,Radius,true);
    auto sc_rc_coca = optimal_reconstruct(dm,vpts_coca,Radius,true);
    auto sc_rc_mult = optimal_reconstruct(dm,vpts_mult,Radius,true);
    
    auto wpk = window_Pk(Radius,0);
    
    auto cc_hl_uni = correlation_coefficients(sc_dm,sc_hl_uni,wpk);
    auto cc_hl = correlation_coefficients(sc_dm,sc_hl,wpk);
    auto cc_mass = correlation_coefficients(sc_dm,sc_rc_mass,wpk);
    auto cc_cntr = correlation_coefficients(sc_dm,sc_rc_cntr,wpk);
    auto cc_envi = correlation_coefficients(sc_dm,sc_rc_envi,wpk);
    auto cc_spin = correlation_coefficients(sc_dm,sc_rc_spin,wpk);
    auto cc_coca = correlation_coefficients(sc_dm,sc_rc_coca,wpk);
    auto cc_mult = correlation_coefficients(sc_dm,sc_rc_mult,wpk);
    
    
    std::cout << "Cross-correlation coefficient:\n";
    std::cout << "[number picture] default   r_n: "              << cc_hl_uni   << ", E= " << sqrt(1-pow(cc_hl_uni,2))  << "\n";
    std::cout << "[mass   picture] default   r_m: "              << cc_hl       << ", E= " << sqrt(1-pow(cc_hl,2))      << "\n";
    std::cout << "---RCST-------Mass:" << Mbin << "-------r_m: " << cc_mass     << ", E= " << sqrt(1-pow(cc_mass,2))    << "\n";
    std::cout << "---RCST-------Envi:" << Ebin << "-------r_m: " << cc_envi     << ", E= " << sqrt(1-pow(cc_envi,2))    << "\n";
    std::cout << "---RCST-------Conc:" << Cbin << "-------r_m: " << cc_cntr     << ", E= " << sqrt(1-pow(cc_cntr,2))    << "\n";
    std::cout << "---RCST-------Dens:" << Dbin << "-------r_m: " << cc_spin     << ", E= " << sqrt(1-pow(cc_spin,2))    << "\n";
    std::cout << "---RCST----M + E + C + D---r_m: "              << cc_coca     << ", E= " << sqrt(1-pow(cc_coca,2))    << "\n";
    std::cout << "---RCST----M * E * C * D---r_m: "              << cc_mult     << ", E= " << sqrt(1-pow(cc_mult,2))    << "\n";
    
    double rmin{5},rmax{150};
    auto vec_r = linear_scale_generator(rmin, rmax, 60, false);
    force_kernel_type(0);
    auto s_rc_mass = sc_back(sc_rc_mass,false);
    auto s_rc_mult = sc_back(sc_rc_mult,false);
    std::vector<double> corr_dm, corr_hl, corr_hl_uni, corr_rc_mass, corr_rc_mult;
    for(int i = 0; auto r : vec_r){
        auto w = wfc(r,0);
        auto c_dm = convol_c2r(sc_dm, w);
        auto c_hl = convol_c2r(sc_hl, w);
        auto c_hl_uni = convol_c2r(sc_hl_uni, w);
        auto c_rc_mass = convol_c2r(sc_rc_mass, w);
        auto c_rc_mult = convol_c2r(sc_rc_mult, w);

        corr_dm.push_back(inner_product(s_dm,c_dm,GridVol)/pow(array_sum(s_dm,GridVol),2)*GridVol -1);
        corr_hl.push_back(inner_product(s_hl,c_hl,GridVol)/pow(array_sum(s_hl,GridVol),2)*GridVol -1);
        corr_hl_uni.push_back(inner_product(s_hl_uni,c_hl_uni,GridVol)/pow(array_sum(s_hl_uni,GridVol),2)*GridVol -1);
        corr_rc_mass.push_back(inner_product(s_rc_mass,c_rc_mass,GridVol)/pow(array_sum(s_rc_mass,GridVol),2)*GridVol -1);
        corr_rc_mult.push_back(inner_product(s_rc_mult,c_rc_mult,GridVol)/pow(array_sum(s_rc_mult,GridVol),2)*GridVol -1);

        delete[] c_dm;
        delete[] c_hl;
        delete[] c_hl_uni;
        delete[] c_rc_mass;
        delete[] c_rc_mult;
        delete[] w;
        std::cout << i << "\n";
        ++i;
    }
    for(auto x : corr_dm) std::cout << x << ", ";std::cout << "\n";
    for(auto x : corr_hl) std::cout << x << ", ";std::cout << "\n";
    for(auto x : corr_hl_uni) std::cout << x << ", ";std::cout << "\n";
    for(auto x : corr_rc_mass) std::cout << x << ", ";std::cout << "\n";
    for(auto x : corr_rc_mult) std::cout << x << ", ";std::cout << "\n";

}


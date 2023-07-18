#include"mracs.h"

void prj_pdf(std::string ofn, double* c, std::vector<Particle>& p0, double min, double max, int nbin);

int main(){
    read_parameter();
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto hl_uni = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    // ----reconstructed hl----
    auto vpts  = halo_mass_split(hl,4);
    auto sc_rc = optimal_reconstruct(dm,vpts,Radius,true);

    // ------convolve with window w------
    auto w = wfc(Radius,0);
    auto s_dm = convol3d(sfc(dm),w,true);
    auto s_hl_uni = convol3d(sfc(hl_uni),w,true);
    auto s_hl = convol3d(sfc(hl),w,true);
    auto s_rc = convol_c2r(sc_rc, w);

    // ----probability distribution function----
    auto p0 = generate_random_particle(500,SimBoxL,0);

    std::vector<std::string> vec_ifn {"dm","hl_uni","hl","rc"};
    std::vector<double*> vec_s {s_dm,s_hl_uni,s_hl,s_rc};
    for(int i = 0; i < vec_ifn.size(); ++i){
        auto ofn = "output/rcpdf_" + vec_ifn[i] + "_" + RADII + ".txt";
        prj_pdf(ofn,vec_s[i],p0,0,5,1000);
    }

}

void prj_pdf(std::string ofn, double* c, std::vector<Particle>& p0, double min, double max, int nbin)
{
    std::ofstream ofs {ofn};
    if(!ofs){
        std::cout << "[pfj_pdf]: opening file " + ofn + " with error, abort!";
        std::terminate();
    }
    const double delta {max - min};
    const double norm {array_sum(c,GridVol) / GridVol};

    auto n = project_value(c,p0,true);

    #pragma omp parallel for
    for(size_t i = 0; i < p0.size(); ++i) n[i] /= norm;

    std::vector<size_t> count(nbin,0);

    //#pragma omp parallel for reduction (+:count)
    for(size_t i = 0; i < p0.size(); ++i) {
        if(min <= n[i] && n[i] < max){
            int index = (n[i] - min)/delta * nbin;
            ++count[index];
        }
    }

    for(int i = 0; i < nbin; ++i){
        ofs << count[i] / delta * nbin / p0.size() << " ";
    }

    ofs << min << " " << max << " ";

}

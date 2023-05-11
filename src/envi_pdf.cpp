// *******************************************************
// print halo environment information at its nearst grid point
// *******************************************************
#include"mracs.h"
using namespace std;

void prj_pdf_temp(std::vector<Particle>& p0, double* c, double nf, double rhomin, double rhomax, int nbin, std::string ofname);


int main(){
    read_parameter();
    std::string ifname {"output/envi_J10_halo_Mcut2e12.txt"};
    std::ifstream ifs {ifname};
    int temp{0};
    char comma{0};
    std::vector<int> envi;
    while(ifs >> temp >> comma) envi.push_back(temp);
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    if(envi.size() == p2.size()) std::cout << "size matched, continue" << std::endl;
    else return 0;

    std::vector<Particle> sheets,filaments,knots;
    for(size_t i = 0; i < p2.size(); ++i)
    {
        if(envi[i] == 1) sheets.push_back(p2[i]);
        else if(envi[i] == 2) filaments.push_back(p2[i]);
        else if(envi[i] == 3) knots.push_back(p2[i]);
    }
    std::cout << "numbers of sheet: " << sheets.size() << std::endl;
    std::cout << "numbers of filament: " << filaments.size() << std::endl;
    std::cout << "numbers of knot: " << knots.size() << std::endl;

    auto p0 = generate_random_particle(1000,SimBoxL,0);

    force_kernel_type(1);
    auto w = wfc(Radius,0);

    auto s = sfc(p1);
    auto c = convol3d(s,w);
    delete[] s;
    prj_pdf_temp(p0, c, static_cast<double>(p1.size())/GridVol,0,5,100,"dm");

    auto s2 = sfc(p2);
    auto c2 = convol3d(s2,w);
    delete[] s2;
    prj_pdf_temp(p0, c2, static_cast<double>(p2.size())/GridVol,0,5,100,"hl");

    auto s_st = sfc(sheets);
    auto c_st = convol3d(s_st,w);
    delete[] s_st;
    prj_pdf_temp(p0, c_st, static_cast<double>(sheets.size())/GridVol,0,5,100,"st");

    auto s_fl = sfc(filaments);
    auto c_fl = convol3d(s_fl,w);
    delete[] s_fl;
    prj_pdf_temp(p0, c_fl, static_cast<double>(filaments.size())/GridVol,0,5,100,"fl");

    auto s_kt = sfc(knots);
    auto c_kt = convol3d(s_kt,w);
    delete[] s_kt;
    prj_pdf_temp(p0, c_kt, static_cast<double>(knots.size())/GridVol,0,5,100,"kt");


}


// nf is the normalization factor, i.e sum(sfc)/(2^J)^3 = expectation of projected value
// e.g. nbin = 100, rhomax = 5, rhomin = 0
void prj_pdf_temp(std::vector<Particle>& p0, double* c, double nf, double rhomin, double rhomax, int nbin, std::string ofname)
{
    std::string ofn = "output/envi_pdf_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofbn = "output/envi_xbin_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::ofstream ofs{ofn}, ofsbin{ofbn};
    if(!ofs || !ofsbin){
        std::cout << "openning file " << ofn << " and output/xbin.txt with error, Abort!" << std::endl;
        std::terminate();
    }
    const double cicexpect = nf;
    const double d_rho = (rhomax - rhomin) / nbin;
    std::cout << "cic expectation: " << cicexpect << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = rhomin + (i + 0.5) * d_rho;
    double count[nbin]{0};
    double value[nbin]{0};

    auto n_prj = project_value(c,p0);
    delete[] c;

    #pragma omp parallel for reduction (+:count)
    for(size_t i = 0; i < p0.size(); ++i) {
        int index = (n_prj[i] / cicexpect - rhomin) / d_rho;
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / p0.size() / d_rho;
    }
    for(int i = 0; i < nbin; ++i) ofsbin << rho[i] << ", "; ofsbin.close();
    for(int i = 0; i < nbin; ++i) ofs << value[i] << ", "; ofs.close();
    std::cout << "========PDF out put========" << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl; 
}
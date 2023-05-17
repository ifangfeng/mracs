// *******************************************************
// print halo environment information at its nearst grid point
// *******************************************************
#include"mracs.h"
using namespace std;

void prj_pdf_temp(std::vector<Particle>& p0, double* c, double R, double nf, double rhomin, double rhomax, int nbin, std::string ofname, bool HM);
std::string GSR {"_GSR5"};

int main(){
    read_parameter();
    std::string ifname {"output/envi_J10" + GSR + "_halo_Mcut2e12.txt"};
    std::ifstream ifs {ifname};
    int temp{0};
    char comma{0};
    std::vector<int> envi;
    while(ifs >> temp >> comma) envi.push_back(temp);
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    if(envi.size() == hl.size()) std::cout << "size matched, continue" << std::endl;
    else {std::cout << "loading data with error, Abort\n"; return 0;}

    std::vector<Particle> sheets,filaments,knots;
    for(size_t i = 0; i < hl.size(); ++i){
        if(envi[i] == 1) sheets.push_back(hl[i]);
        else if(envi[i] == 2) filaments.push_back(hl[i]);
        else if(envi[i] == 3) knots.push_back(hl[i]);
    }

    std::vector<double> vecR {20,30,40,50,80};
    std::vector<std::vector<Particle>> vecP {dm, hl, sheets, filaments, knots};
    std::vector<string> ofname {"dm", "hl", "st", "fl", "kt"};
    for(int i = 0; i < vecP.size(); ++i){
        std::cout << "number of " << ofname[i] << ": " << vecP[i].size() << "\n";
    }

    force_kernel_type(1);
    auto p0 = generate_random_particle(1000,SimBoxL,0);
    
    for(auto R : vecR){
        auto w = wfc(R,0);
        for(int i = 0; i < vecP.size(); ++i){
            auto s = sfc(vecP[i]);
            auto c = convol3d(s,w);
            double norm = static_cast<double>(vecP[i].size())/GridVol;
            delete[] s;
            int xbin {200};
            if(R < 40 && i == 4) xbin = 100;
            if(R < 30 && (i != 0)) xbin = 50; 
            double rhomax = 5;
            if(R > 30) rhomax = 3;
            prj_pdf_temp(p0,c,R,norm,0,rhomax,xbin,ofname[i],1);
        }
    }

}


// nf is the normalization factor, i.e sum(sfc)/(2^J)^3 = expectation of projected value
// e.g. nbin = 100, rhomax = 5, rhomin = 0
void prj_pdf_temp(std::vector<Particle>& p0, double* c, double R, double nf, double rhomin, double rhomax, int nbin, std::string ofname, bool HM)
{
    std::string suffix = ofname + "_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution);
    if(HM) suffix += "_HM" + GSR + ".txt";
    else suffix += GSR + ".txt";
    std::string ofn = "output/envi_pdf_" + suffix;
    std::string ofbn = "output/envi_xbin_" +suffix;
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
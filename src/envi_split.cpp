// *******************************************************
// halo environment split and sample for scatter plot
// *******************************************************
#include"mracs.h"

std::string GSR {"_GSR5"};

int main(){
    read_parameter();
    std::string ifname {"output/envi_J10" + GSR + "_halo_Mcut2e12.txt"};
    std::ifstream ifs {ifname};
    std::vector<int> envi;int temp{0}; char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);

    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    if(envi.size() == hl.size()) std::cout << "size matched, continue" << std::endl;
    else {std::cout << "loading data with error, Abort\n"; return 0;}

    std::vector<Particle> voids,sheets,filaments,knots;
    for(size_t i = 0; i < hl.size(); ++i){
        if(envi[i] == 0) voids.push_back(hl[i]);
        else if(envi[i] == 1) sheets.push_back(hl[i]);
        else if(envi[i] == 2) filaments.push_back(hl[i]);
        else if(envi[i] == 3) knots.push_back(hl[i]);
    }

    auto p0 = generate_random_particle(20,SimBoxL,0);
    std::vector<double> vecR {20,30,40,50,80};
    std::vector<std::vector<Particle>> vecP {dm, hl, voids, sheets, filaments, knots};
    std::vector<double*> svec;
    for(auto x : vecP) svec.push_back(sfc(x));

    std::string prefix {"output/envi_split_scatter_"};
    std::vector<std::string> fname {"dm","hl","vd","st","fl","kt"};
    std::string suffix {"_HM_halo_Mcut2e12.txt"};
    
    force_kernel_type(1);
    for(auto R : vecR){
        auto w = wfc(R,0);
        for(int i = 0; i < vecP.size(); ++i){
            auto c = convol3d(svec[i],w);
            auto n = project_value(c,p0,true);
            double cicexpect = static_cast<double>(vecP[i].size())/GridVol;
            std::string ofname = prefix + fname[i] +"_R" + std::to_string((int)R) + "_J" + std::to_string(Resolution) + suffix;
            std::ofstream ofs {ofname};
            for(size_t j = 0; j < p0.size(); ++j) ofs << n[j]/cicexpect << ", ";
            delete[] n;
        }
    }
}



// c is cic counting result, n is prj result, rhomin and rhomax could be 0, 5 respectively
void cp_temp_dispersion(std::vector<int64_t>& c, double* n, double rhomin, double rhomax, double cicexpect, std::string ofname)
{
    std::string suffix = "_R" +  std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofxbin = "output/cp_xbin_" + ofname + suffix;
    std::string ofmean = "output/cp_mean_" + ofname + suffix;
    std::string ofrdev = "output/cp_rdev_" + ofname + suffix;
    std::string ofrmean = "output/cp_rmean_" + ofname + suffix;
    std::string ofdev = "output/cp_stddev_" + ofname + suffix;
    std::string ofxs = "output/cp_scatter_x_" + ofname + suffix;
    std::string ofys = "output/cp_scatter_y_" + ofname + suffix;
    std::string ofrys = "output/cp_scatter_ry_" + ofname + suffix;

    std::ofstream ofsxbin{ofxbin}, ofsmean{ofmean}, ofsdev{ofdev}, ofsrmean{ofrmean}, ofsrdev{ofrdev}, ofsxs{ofxs}, ofsys{ofys}, ofsrys{ofrys};
    if(!ofsxbin || !ofsmean || !ofsdev ||!ofsrmean || !ofsrdev || !ofsxs || !ofsys || !ofsrys){
        std::cout << "openning file " << ofmean << " and " << ofdev << " with error, Abort!" << std::endl;
        std::terminate();
    }
    const int nscatter = 10000;
    const int increment = c.size() / nscatter;
    const int npt = (rhomax - rhomin) * cicexpect + 1;
    const int xresolmax = 20;
    const int xstep = cicexpect / xresolmax + 1;
    const int nbin = npt / xstep;
    double rhox[npt]; for(int i = 0; i < npt; ++i) rhox[i] = i / cicexpect;
    double xbin[nbin]; for(int i = 0; i < nbin; ++i) xbin[i] = i / cicexpect * xstep; // delta_xbin = xstep / cicexpect
    double count[nbin]{0};
    double mean[nbin]{0};
    double rmean[nbin]{0};
    double deviation[nbin]{0};
    double relatdev[nbin]{0};

    #pragma omp parallel for reduction (+:mean,deviation,count)
    for(size_t i = 0; i < c.size(); ++i){
        if(c[i] < npt){
            mean[c[i]/xstep] += n[i];
            deviation[c[i]/xstep] += n[i] * n[i];
            ++count[c[i]/xstep];
        }
    }
    for(int i = 0; i < nbin; ++i){
        if(count[i] > 1){
            mean[i] /= count[i];
            deviation[i] = sqrt((deviation[i] - pow(mean[i],2) * count[i]) / (count[i] - 1));
            relatdev[i] = (i ? (deviation[i] / (i*xstep+xstep/2)):1);
            rmean[i] = (i ? ((mean[i] - (i*xstep+xstep/2)) / (i*xstep+xstep/2)):1);
        }
        else if(count[i] == 1){
            rmean[i] = (i ? ((mean[i] - (i*xstep+xstep/2)) / (i*xstep+xstep/2)):1);
        }
    }
    for(int i = 0; i < nbin; ++i) ofsxbin << xbin[i] << ", "; ofsxbin.close();
    for(int i = 0; i < nbin; ++i) ofsmean << (mean[i] ? (mean[i] - (i*xstep+xstep/2)) : 0) << ", "; ofsmean.close();
    for(int i = 0; i < nbin; ++i) ofsrmean << ((rmean[i] < 1)? rmean[i] : 1) << ", "; ofsmean.close();
    for(int i = 0; i < nbin; ++i) ofsdev << deviation[i] << ", "; ofsdev.close();
    for(int i = 0; i < nbin; ++i) ofsrdev << ((relatdev[i] < 1)? relatdev[i] : 1) << ", "; ofsrdev.close();
    for(int i = 0; i < c.size(); i += increment){
        if(c[i] < npt){
            ofsxs << c[i] / cicexpect << ", ";
            ofsys << n[i] - c[i] << ", ";
            ofsrys << (c[i]? (((n[i] - c[i]) / c[i] > 1)? 1 : ((n[i] - c[i]) / c[i])) : 1) << ", ";
        }
    }
    ofsxs.close();
    ofsys.close();
}
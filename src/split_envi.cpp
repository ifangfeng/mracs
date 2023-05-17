// *******************************************************
// halo environment split and sample for scatter plot
// *******************************************************
#include"mracs.h"

int main(){
    read_parameter();

    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto p0 = generate_random_particle(100,SimBoxL,0);

    double GSR {5};
    force_base_type(0,1);
    force_kernel_type(2);
    auto s = sfc(p1);
    auto sc = sfc_r2c(s);
    auto w = wft(GSR, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p0);
    delete[] s;
    delete[] w;
    delete[] sc;

    std::vector<double> vecR {20,30,40,50,80};
    force_base_type(0,5);
    force_kernel_type(1);
    auto s1 = sfc(p1);
    auto s2 = sfc(p2);
    auto sc1 = sfc_r2c(s1);
    auto sc2 = sfc_r2c(s2);

    std::vector<double> voids_x, sheets_x, filaments_x, knots_x;
    std::vector<double> voids_y, sheets_y, filaments_y, knots_y;

    for(auto R : vecR){
        auto w = wfc(R,0);
        auto c1 = convol_c2r(sc1,w);
        auto c2 = convol_c2r(sc2,w);
        auto n1 = project_value(c1,p0);
        auto n2 = project_value(c2,p0);

        for(int i = 0; i < env.size(); ++i){
            if (env[i] == 0) {
                voids_x.push_back(n1[i]);
                voids_y.push_back(n2[i]);}
            else if (env[i] == 1){
                sheets_x.push_back(n1[i]);
                sheets_y.push_back(n2[i]);}
            else if (env[i] == 2){
                filaments_x.push_back(n1[i]);
                filaments_y.push_back(n2[i]);}
            else {
                knots_x.push_back(n1[i]);
                knots_y.push_back(n2[i]);}
        }

    }


    


    int64_t voids{0}, sheets{0}, filaments{0}, knots{0};
    for(auto x : env){
        
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
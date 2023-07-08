// =======================================================
// This cpp source file is part of MRACS project
// one point statistics
// =======================================================
#include"opstat.h"



// nf is the normalization factor, i.e sum(sfc)/(2^J)^3 = expectation of projected value
// e.g. nbin = 100, rhomax = 5, rhomin = 0
void pdf(std::vector<Particle>& p0, double* c, double nf, double rhomin, double rhomax, int nbin, std::string ofname)
{
    std::string ofn = "output/prj_pdf_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofn1 = "output/prj_M1th_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofn2 = "output/prj_M2nd_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofn3 = "output/prj_M3rd_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::string ofbn = "output/prj_xbin_" + ofname + "_R" + std::to_string((int)Radius) + "_J" + std::to_string(Resolution) + ".txt";
    std::ofstream ofs{ofn}, ofsbin{ofbn}, ofs1{ofn1}, ofs2{ofn2}, ofs3{ofn3};
    if(!ofs || !ofsbin || !ofs1 || !ofs2 || !ofs3){
        std::cout << "openning file " << ofn << " and output/xbin.txt with error, Abort!" << std::endl;
        std::terminate();
    }
    const double cicexpect = nf;
    const double d_rho = (rhomax - rhomin) / nbin;
    std::cout << "cic expectation: " << cicexpect << std::endl;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = rhomin + (i + 0.5) * d_rho;
    double count[nbin]{0};
    double value[nbin]{0};
    double mmt1[nbin]{0};
    double mmt2[nbin]{0};
    double mmt3[nbin]{0};
    auto n_prj = project_value(c,p0,true);

    #pragma omp parallel for reduction (+:count)
    for(size_t i = 0; i < p0.size(); ++i) {
        int index = (n_prj[i] / cicexpect - rhomin) / d_rho;
        if(index < nbin && index >= 0)
        ++count[index];
    }
    for(int i = 0; i < nbin; ++i){
        value[i] = count[i] / p0.size() / d_rho;
        mmt1[i] = rho[i] * value[i];
        mmt2[i] = rho[i] * rho[i] * value[i];
        mmt3[i] = rho[i] * rho[i] * rho[i] * value[i];
    }
    for(int i = 0; i < nbin; ++i) ofsbin << rho[i] << ", "; ofsbin.close();
    for(int i = 0; i < nbin; ++i) ofs << value[i] << ", "; ofs.close();
    for(int i = 0; i < nbin; ++i) ofs1 << mmt1[i] << ", "; ofs1.close();
    for(int i = 0; i < nbin; ++i) ofs2 << mmt2[i] << ", "; ofs2.close();
    for(int i = 0; i < nbin; ++i) ofs3 << mmt3[i] << ", "; ofs3.close();
    std::cout << "========PDF out put========" << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl; 
}

// c is cic counting result, n is prj result, rhomin and rhomax could be 0, 5 respectively
void cp_dispersion(std::vector<int64_t>& c, double* n, double rhomin, double rhomax, double cicexpect, std::string ofname)
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

    //std::cout << "========cic PDF out put========" << std::endl;
    //std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    //std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    //std::cout << "mean: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl; 
    //std::cout << "variance: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << variance[i] << ", "; std::cout << std::endl; 
    //std::cout << "rv: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << relatdev[i] << ", "; std::cout << std::endl; 
}


// c is CIC sampling result , other value would be looks like: rhomax = 5, rhomin = 0
void cic_pdf(std::vector<int64_t>& c, double rhomin, double rhomax, double cicexpect, std::string ofname)
{
    std::string ofn = "output/cic_pdf_" + ofname + "_R" + std::to_string((int)Radius) + ".txt";
    std::string ofbn = "output/cic_xbin_" + ofname + "_R" +  std::to_string((int)Radius) + ".txt";
    std::string ofn1 = "output/cic_M1th_" + ofname + "_R" + std::to_string((int)Radius) + ".txt";
    std::string ofn2 = "output/cic_M2nd_" + ofname + "_R" + std::to_string((int)Radius) + ".txt";
    std::string ofn3 = "output/cic_M3rd_" + ofname + "_R" + std::to_string((int)Radius) + ".txt";
    std::ofstream ofs{ofn}, ofsbin{ofbn}, ofs1{ofn1}, ofs2{ofn2}, ofs3{ofn3};
    if(!ofs || !ofsbin || !ofs1 || !ofs2 || !ofs3){
        std::cout << "openning file " << ofn << " and " << ofbn << " with error, Abort!" << std::endl;
        std::terminate();
    }
    const int nbin = (rhomax - rhomin) * cicexpect + 1;
    double rho[nbin]; for(int i = 0; i < nbin; ++i) rho[i] = i / cicexpect;
    double count[nbin]{0};
    double value[nbin]{0};
    double mmt1[nbin]{0};
    double mmt2[nbin]{0};
    double mmt3[nbin]{0};

    #pragma omp parallel for reduction (+:count)
    for(size_t i = 0; i < c.size(); ++i){
        if(c[i] < nbin)
            ++count[c[i]];
    }
    for(size_t i = 0; i < nbin; ++i){
        value[i] = count[i] / c.size() * cicexpect; // delta_rho = 1. / cicecpect
        mmt1[i] = i * count[i] / c.size();
        mmt2[i] = i * rho[i] * count[i] / c.size();
        mmt3[i] = i * rho[i] * rho[i] * count[i] / c.size();
    }
    for(int i = 0; i < nbin; ++i) ofsbin << rho[i] << ", "; ofsbin.close();
    for(int i = 0; i < nbin; ++i) ofs << value[i] << ", "; ofs.close();
    for(int i = 0; i < nbin; ++i) ofs1 << mmt1[i] << ", "; ofs1.close();
    for(int i = 0; i < nbin; ++i) ofs2 << mmt2[i] << ", "; ofs2.close();
    for(int i = 0; i < nbin; ++i) ofs3 << mmt3[i] << ", "; ofs3.close();

    std::cout << "========cic PDF out put========" << std::endl;
    std::cout << "density bin: " << std::endl;  for(int i = 0; i < nbin; ++i) std::cout << rho[i] << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl;        for(int i = 0; i < nbin; ++i) std::cout << count[i] << ", "; std::cout << std::endl;
    std::cout << "profile: " << std::endl;      for(int i = 0; i < nbin; ++i) std::cout << value[i] << ", "; std::cout << std::endl; 

}
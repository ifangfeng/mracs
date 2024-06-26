// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

double* fourier_mode_correlation
(std::vector<std::vector<Particle>>& dm,std::vector<std::vector<Particle>>& hl);

int main(){
    read_parameter();
    
    std::vector<double> we{-0.973741, 0.00269408, 0.136522, 0.182164}; //weight
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto halo = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");
    std::string ifname {"output/envi_J10_GSR5_halo_Mcut2e12.txt"};

    std::ifstream ifs {ifname};
    std::vector<int> envi;int temp{0}; char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(halo.size() == envi.size()) std::cout << "halo size matched! continue\n"; else std::terminate();
    std::vector<double> npartenvi(4,0);
    for(auto x : envi) npartenvi[x]++;
    std::vector<Particle> p2;
    for(size_t i = 0; i < envi.size(); ++i){
        p2.push_back({halo[i].x,halo[i].y,halo[i].z,we[envi[i]]*npartenvi[envi[i]]});
    }

    int Nl {4};
    int Nrl {Nl * Nl * Nl}; // number of realization

    std::vector<std::vector<Particle>> dm(Nrl), hl(Nrl);
    SimBoxL /= Nl;
    for(auto x : p1) {
        int i,j,k;
        i = x.x/SimBoxL;
        j = x.y/SimBoxL;
        k = x.z/SimBoxL;
        dm[i * Nl * Nl + j * Nl + k]
        .push_back({x.x - i * SimBoxL, x.y - j * SimBoxL, x.z - k * SimBoxL, 1.});
    }
    for(auto x : p2) {
        int i,j,k;
        i = x.x/SimBoxL;
        j = x.y/SimBoxL;
        k = x.z/SimBoxL;
        hl[i * Nl * Nl + j * Nl + k]
        .push_back({x.x - i * SimBoxL, x.y - j * SimBoxL, x.z - k * SimBoxL, x.weight});
    }

    //std::ofstream ofs0 {"output/data0.txt"};
    std::ofstream ofs1 {"output/data1.txt"};

    //force_base_type(0,1);
    //auto r0 = fourier_mode_correlation(dm, hl);
    //for(size_t i = 0; i < (GridLen/2 + 1) * GridLen * GridLen; ++i) ofs0 << r0[i] << " ";delete[] r0;
    
    force_base_type(1,10);
    auto r1 = fourier_mode_correlation(dm, hl);
    for(size_t i = 0; i < (GridLen/2 + 1) * GridLen * GridLen; ++i) ofs1 << r1[i] << " ";delete[] r1;
    
    
}

// r(k), k is an vector. correlation as function of fourier mode of two density fields dm and hl
// Notice
double* fourier_mode_correlation
(std::vector<std::vector<Particle>>& dm,std::vector<std::vector<Particle>>& hl)
{
    double* mm = new double[(GridLen/2 + 1) * GridLen * GridLen];
    double* hh = new double[(GridLen/2 + 1) * GridLen * GridLen];
    double* mh = new double[(GridLen/2 + 1) * GridLen * GridLen];

    for(int i = 0; i < dm.size(); ++i){
        auto dm_sc = sfc_r2c(sfc(dm[i]),true);
        auto hl_sc = sfc_r2c(sfc(hl[i]),true);

        #pragma omp parallel for
        for(size_t l = 0; l < (GridLen/2 + 1) * GridLen * GridLen; ++l){
            mm[l] += pow(dm_sc[l][0],2) + pow(dm_sc[l][1],2);
        }
        #pragma omp parallel for
        for(size_t l = 0; l < (GridLen/2 + 1) * GridLen * GridLen; ++l){
            hh[l] += pow(hl_sc[l][0],2) + pow(hl_sc[l][1],2);
        }
        #pragma omp parallel for
        for(size_t l = 0; l < (GridLen/2 + 1) * GridLen * GridLen; ++l){
            mh[l] += dm_sc[l][0] * hl_sc[l][0] + dm_sc[l][1] * hl_sc[l][1];
        }
        fftw_free(dm_sc);
        fftw_free(hl_sc);
    }

    #pragma omp parallel for
    for(size_t l = 0; l < (GridLen/2 + 1) * GridLen * GridLen; ++l){
        mh[l] /= sqrt(mm[l] * hh[l]);
    }
    delete[] mm;
    delete[] hh;

    return mh;
}

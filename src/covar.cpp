// covariance of different halo bin 
#include"mracs.h"

void halo_envi_match(std::string ifn, std::vector<Particle>& hl, std::vector<Particle>& vd
, std::vector<Particle>& st, std::vector<Particle>& fl, std::vector<Particle>& kt);


int main(){
    read_parameter();
    
    auto dms = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hls = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    std::string ifname {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    std::vector<Particle> vds,sts,fls,kts;
    halo_envi_match(ifname,hls,vds,sts,fls,kts);
    std::vector<std::vector<Particle>*> vpts {&dms,&hls,&vds,&sts,&fls,&kts};

    std::vector<fftw_complex*> vec_sc;
    for(auto x : vpts) vec_sc.push_back(sfc_r2c(sfc(*x),true));

    double rmin{100},rmax{100}; 
    auto vecR = log_scale_generator(rmin,rmax,1,false);
    std::vector<double*> vec_wpk; for(auto r : vecR) vec_wpk.push_back(window_Pk(r,0));
    
    std::vector<double> cov,corr;
    std::vector<double> weight;
    for(int i = 2; i < vpts.size(); ++i) weight.push_back((*vpts[i]).size()/static_cast<double>((*vpts[1]).size()));
    double deno{0};for(auto x : weight) deno+=x;double we[4];for(int i = 0; i < 4; ++i) we[i] = weight[i]; 
    for(auto x : weight) std::cout << x << ", "; std::cout << deno << std::endl;


    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = 0; j <= i; ++j)
            for(auto w : vec_wpk){
                cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),w));
            }
    std::cout << "covariance: \n";
    std::cout << "{array element: (dm,hl,vd,st,fl,kt)^T x (dm,hl,vd,st,fl,kt)}\n";
    for(int i = 0; i < vec_sc.size(); ++i){
        for(int j = 0; j <= i; ++j)
            std::cout << cov[i*(i+1)/2 + j] << ", ";
        std::cout << std::endl;
    }
    double tmp{0};
    for(int i = 2; i < vec_sc.size(); ++i)
        tmp += weight[i-2] * cov[i*(i+1)/2];
    std::cout << "predict: " << tmp << ", calculate: " << cov[1] << std::endl;


    std::cout << "cross-correlation: \n";
    std::cout << "(Array element: {(dm,hl,vd,st,fl,kt)^T x (dm,hl,vd,st,fl,kt)})\n";
    for(int i = 0; i < vec_sc.size(); ++i){
        for(int j = 0; j <= i; ++j)
            std::cout << cov[i*(i+1)/2 + j]/sqrt(cov[i*(i+1)/2 + i]*cov[j*(j+1)/2 + j]) << ", ";
        std::cout << std::endl;
    }
    
    auto array = new double[16];
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j)
        {
            if(j > i){
                array[i * 4 + j] = cov[(j+2)*(j+3)/2 + i+2];
            }
            else 
            {
                array[i * 4 + j] = cov[(i+2)*(i+3)/2 + j+2];
            }
        }
    }
    //std::cout << "matrix A: \n";
    //for(int i = 0; i < 4; ++i){
    //    for(int j = 0; j < 4; ++j){
    //        std::cout << array[i * 4 + j] << ", ";
    //    }std::cout << std::endl;
    //}
    double array_b[4];
    for(int i = 0; i < 4; ++i) array_b[i] = cov[(i+2)*(i+3)/2];
    //std::cout << "vector beta: \n";
    //for(int i = 0; i < 4; ++i) std::cout << array_b[i] << ", "; std::cout << std::endl;

    std::cout << "default:\n";
    std::cout << "w: " << we[0] << ", "  << we[1] << ", " << we[2] << ", " << we[3] << std::endl;
    double ccr{cov[1]/sqrt(cov[2])}; std::cout << "r: " << ccr / sqrt(cov[0]) << std::endl;
    int N{1000};
    int I{0}, J{0}, K{0};

    double sum{0};
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
            {
                double w[4];
                w[0] = (N - i - j - k)/static_cast<double>(N);
                w[1] = i/static_cast<double>(N);
                w[2] = j/static_cast<double>(N);
                w[3] = k/static_cast<double>(N);

                sum = 0;
                for(int m = 0; m < 4; ++m)
                    for(int n = 0; n < 4; ++n)
                        sum += array[m * 4 + n] * w[m] * w[n];
                double r = (array_b[0] * w[0] + array_b[1] * w[1] + array_b[2] * w[2] + array_b[3] * w[3]) / sqrt(sum);
                if (r > ccr) {ccr = r; I = i; J = j; K = k;}
            }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Param Explore T = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
    << "[ms]" << std::endl;

    std::cout << "Max: \n";

    double w[4];
    w[0] = (N - I - J - K)/static_cast<double>(N);
    w[1] = I/static_cast<double>(N);
    w[2] = J/static_cast<double>(N);
    w[3] = K/static_cast<double>(N);
    double sum0{0};
    for(int m = 0; m < 4; ++m)
        for(int n = 0; n < 4; ++n)
            sum0 += array[m * 4 + n] * w[m] * w[n];
    double r = (array_b[0] * w[0] + array_b[1] * w[1] + array_b[2] * w[2] + array_b[3] * w[3]) / sqrt(sum0 * cov[0]);

    std::cout << "w: " << w[0] << ", " << w[1] << ", " << w[2] << ", " << w[3] << std::endl;
    std::cout << "r: " << r << std::endl;

    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));

    for(size_t n = 0; n < GridLen * GridLen * (GridLen/2 + 1); ++n){
        sc[n][0] = vec_sc[2][n][0] * w[0] / we[0] + vec_sc[3][n][0] * w[1] / we[1] +
                   vec_sc[4][n][0] * w[2] / we[2] + vec_sc[5][n][0] * w[3] / we[3];
        sc[n][1] = vec_sc[2][n][1] * w[0] / we[0] + vec_sc[3][n][1] * w[1] / we[1] +
                   vec_sc[4][n][1] * w[2] / we[2] + vec_sc[5][n][1] * w[3] / we[3];
    }
    double ab = covar_CombinewithKernel(densityCovarianceArray(sc,vec_sc[0]),vec_wpk[0]);
    double aa = covar_CombinewithKernel(densityCovarianceArray(sc,sc),vec_wpk[0]);

    std::cout << "r(reconstruct): " << ab/sqrt(aa * cov[0]) << std::endl;

/*
    auto p0 = default_random_particle(SimBoxL,1000*10);
    auto win = wfc(Radius,0);
    auto n = project_value(convol_c2r(sc,win),p0,true);
    auto n0= project_value(convol_c2r(vec_sc[1],win),p0,true);
    auto nm= project_value(convol_c2r(vec_sc[0],win),p0,true);

    //std::cout << sc[0][0] << " and " << hls.size() << std::endl;

    std::ofstream ofs {"output/scatter_cp_"+RADII+".txt"};
    for(size_t i = 0; i < p0.size(); ++i) 
    ofs << nm[i]/dms.size()*GridVol << " " << n0[i]/hls.size()*GridVol << " " << n[i]/sc[0][0]*GridVol << " ";
    
    double wpca[4]{-8.81084, 0.0134186, 0.751282, 1};
*/
}



void halo_envi_match(std::string ifn, std::vector<Particle>& hl, std::vector<Particle>& vd
, std::vector<Particle>& st, std::vector<Particle>& fl, std::vector<Particle>& kt)
{
    std::ifstream ifs {ifn};
    if(!ifs){std::cout << "reading " + ifn + " with error, Abort"; std::terminate();}

    std::vector<int> envi;
    int temp{0};
    char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) 
        std::cout << "halo size matched! continue\n"; 
    else std::terminate();

    for(size_t i = 0; i < envi.size(); ++i){
        if(envi[i] == 0) vd.push_back(hl[i]);
        else if(envi[i] == 1) st.push_back(hl[i]);
        else if(envi[i] == 2) fl.push_back(hl[i]);
        else if(envi[i] == 3) kt.push_back(hl[i]);
    }
    std::vector<int>().swap(envi);
    if(hl.size() != (vd.size() + st.size() + fl.size() + kt.size())) 
        std::cout << "Warning! halo environment subset size not matched\n";
}
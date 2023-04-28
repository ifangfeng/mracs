#include"mracs.h"

double* prj_test1(const double* s, std::vector<Particle>& p0);
double* prj_test2(const double* s, std::vector<Particle>& p0);


int main(){
    read_parameter();
    std::vector<Particle> p = read_in_DM_3vector(DataDirec);
    auto s  = sfc(p);
    auto w  = wfc(Radius,0);
    auto c  = convol3d(s,w);

    auto x = prj_test1(c,p);
    auto x2 = prj_test2(c,p);

    for(int i = 0; i < 100; ++i) std::cout << x[i] << "  vs. " << x2[i] << std::endl; 

}



double* prj_test1(const double* s, std::vector<Particle>& p0)
{
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates
    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

    auto begin4 = std::chrono::steady_clock::now();

    

    auto result = new double[p0.size()];
    double sum {0};
    #pragma omp parallel for reduction (+:sum)
    for(size_t n = 0; n < p0.size(); ++n)
    {
        double xx = p0[n].x * ScaleFactor;
        double yy = p0[n].y * ScaleFactor;
        double zz = p0[n].z * ScaleFactor;

        int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
        int zzc = floor(zz), zzf = SampRate * (zz - zzc);

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                {
                    sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                            * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }

        result[n] = sum;
        sum = 0;
    }

    auto end4 = std::chrono::steady_clock::now();
    std::cout << "Test1 = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return result;
}

double* prj_test2(const double* s, std::vector<Particle>& p0)
{
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates

    auto begin4 = std::chrono::steady_clock::now();

    auto result = new double[p0.size()];
    double sum {0};
    #pragma omp parallel for reduction (+:sum)
    for(size_t n = 0; n < p0.size(); ++n)
    {
        double xx = p0[n].x * ScaleFactor;
        double yy = p0[n].y * ScaleFactor;
        double zz = p0[n].z * ScaleFactor;

        int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
        int zzc = floor(zz), zzf = SampRate * (zz - zzc);

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                {
                    sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                            * phi[xxf + i*SampRate] * phi[yyf + j*SampRate] * phi[zzf + k*SampRate];
                }

        result[n] = sum;
        sum = 0;
    }

    auto end4 = std::chrono::steady_clock::now();
    std::cout << "Test2 = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return result;
}
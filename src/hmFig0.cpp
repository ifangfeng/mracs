#include"mracs.h"

double* sfc_tmp2d(std::vector<Particle>& p);
std::vector<double> proj_tmp2d(const double* s, std::vector<Particle>& p0);

int main(){
    read_parameter();
    //auto p = default_random_particle(SimBoxL,5);
    std::vector<Particle> p; p.push_back({20,20,0,1});p.push_back({26,26,0,1});p.push_back({50,72,0,1});
    p.push_back({68,32,0,1});p.push_back({88,96,0,1});

    std::vector<Particle> p0;
    const int ng{800};
    for(int i = 0; i < ng; ++i){
        for(int j = 0; j < ng; ++j){
            p0.push_back({SimBoxL/ng * i, SimBoxL/ng * j, 0., 1});
        }
    }
    //for(auto x : p0) std::cout << x.x << ". " << x.y << ", ";
    //for(auto x : p) std::cout << x.x << ", " << x.y << ", " << x.z << std::endl;
    std::vector<int> resol{2,3,4,5,6,7};
    for(auto J : resol){
        std::string ofn {"output/hmFig0J"+std::to_string(J)+".txt"};
        std::ofstream ofs{ofn};
        force_resoluton_J(J);
        auto s = sfc_tmp2d(p);
        auto rho = proj_tmp2d(s,p0);
        for(int i = 0; i < p0.size(); ++i) ofs << pow(2,2*J-4)*rho[i] << " ";
    }
}

//=======================================================================================
//||||||||||||||| sfc3d (scaling function coefficients of density field) ||||||||||||||||
//=======================================================================================
double* sfc_tmp2d(std::vector<Particle>& p)
{   
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates

    auto s = new double[GridLen * GridLen]();         // density field coefficients in v_j space

    #pragma omp parallel for
    for (size_t n = 0; n < p.size(); ++n) {
        // rescale the particle coordinates to MRA framework
        double xx = (p[n].x ) * ScaleFactor;
        double yy = (p[n].y ) * ScaleFactor;

        // get the 'Coarse' and 'Finer' coodinate, e.g. if
        // position p == 63.25678, then c == 63, f == 256
        int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy), yyf = SampRate * (yy - yyc);      

        double s_temp[phiSupport * phiSupport]{0};
        
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                s_temp[i*phiSupport + j] += phi[xxf + i * SampRate] * phi[yyf + j * SampRate] * p[n].weight;// / total;

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                s[((xxc - i)&(GridLen - 1))*GridLen + ((yyc - j)&(GridLen - 1))] += s_temp[i*phiSupport + j];
    }
    return s;
}

// given scaling function coefficients s, return the value it represented at position p0
std::vector<double> proj_tmp2d(const double* s, std::vector<Particle>& p0)
{
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates

    std::vector<double> result(p0.size());
    double sum {0};
    #pragma omp parallel for reduction (+:sum)
    for(size_t n = 0; n < p0.size(); ++n)
    {
        double xx = p0[n].x * ScaleFactor;
        double yy = p0[n].y * ScaleFactor;

        int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy), yyf = SampRate * (yy - yyc);      

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
            {
                sum += s[((xxc-i) & (GridLen-1)) * GridLen + ((yyc-j) & (GridLen-1))]
                        * phi[xxf + i * SampRate] * phi[yyf + j * SampRate];
            }

        result[n] = sum;
        sum = 0;
    }

    return result;
}
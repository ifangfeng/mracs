#include"mracs.h"

double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport);
double Proj_Value2(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport);
void Proj_Value3( std::vector<Particle>& p0, double* s, std::vector<int> step, int phiSupport,std::vector<double>& result);
void result_interpret1(const double* s, std::vector<Particle>& p0, std::vector<double>& result);


//--------------------------------------------------------
int main(){

    read_parameter();
    const int len = 1e6;

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0,GridLen);

    std::chrono::steady_clock::time_point generate = std::chrono::steady_clock::now();

    std::vector<Particle> p0,p1;

    for(size_t i = 0; i < len ; ++i) p0.push_back({v(e), v(e), v(e), 1.});
    for(int i = 0; i < 100; i++)
        for(int j = 0; j < 100; j++)
            for(int k = 0; k < 100; k++)
                p1.push_back({(i+u(e))/100*GridLen, (j+u(e))/100*GridLen, (k+u(e))/100*GridLen, 1.});

    double mean {0};
    for(Particle i : p0) mean += i.x;
    std::cout<< mean << std::endl;


    std::chrono::steady_clock::time_point generate2 = std::chrono::steady_clock::now();
    std::cout << "Time difference genaration    = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(generate2 - generate ).count()
    << "[ms]" << std::endl;

    auto s0 = sfc(p0);
    auto w  = wfc(Radius,0);
    auto c0 = convol3d(s0, w);


    const int phiSupport = phi[phi.size() - 1] - phi[phi.size() - 2];
    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i) step[i] = i * SampRate;

    std::cout<< '\n' << std::endl;

//111111111111------------------------------------------------------

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    for (Particle i : p0){
        auto x = Proj_Value2(i.x,i.y,i.z, c0, step, phiSupport);
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "1:Time difference (i:p0)  Proj_Value2    = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1 ).count()
    << "[ms]" << std::endl;


//22222222222222222------------------------------------------------------
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto ptest  = Particle{v(e), v(e), v(e), 1.};
    for (int i = 0 ; i < 1e6 ; i++)
        Proj_Value2(ptest.x,ptest.y,ptest.z,c0,step,phiSupport);


    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "2:Time difference (outside loop)   Proj_Value2   = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2 ).count()
    << "[ms]" << std::endl;

//333333333333333333------------------------------------------------------

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();

    Proj_Value(ptest.x,ptest.y,ptest.z,c0,step,phiSupport);


    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    std::cout << "3:Time difference (inside loop)    Proj_Value   = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3 ).count()
    << "[ms]" << std::endl;

// //4444444444444444------------------------------------------------------




    // std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();
    // std::cout<< len <<std::endl;

    // std::vector<double> projNum3(len);
    // for (int i =0;i<len;i++) projNum3[i] = Proj_Value2(p0[i].x,p0[i].y,p0[i].z, c0, step, phiSupport);



    // std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    // std::cout << "4:Time difference (p0 inside)  projNum_vetcor   = "
    // << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4 ).count()
    // << "[ms]" << std::endl;

//5555555555555555------------------------------------------------------
    
    std::chrono::steady_clock::time_point begin5 = std::chrono::steady_clock::now();

    std::cout<< len <<std::endl;
    double projNum6[len];
    // double* projNum7 = new double [len];
    for (int i =0; i<len;i++) {
        projNum6[i] = Proj_Value2(p0[i].x,p0[i].y,p0[i].z, c0, step, phiSupport);  
    }

    std::chrono::steady_clock::time_point end5 = std::chrono::steady_clock::now();
    std::cout << "5:Time difference (p0 inside)  projNum_array   = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - begin5 ).count()
    << "[ms]" << std::endl;
    
    int ii {0};
    for(int i=0;i<len;i++) {
        if (projNum6[i]<0) ii++;
    } 
    std::cout<< ii <<std::endl;

//6666666666666666------------------------------------------------------
    // std::vector<double> projNum4;
    // result_interpret1(c0, p0, projNum4);



//666666666666666------------------------------------------------------
    
    std::chrono::steady_clock::time_point begin6 = std::chrono::steady_clock::now();

    std::cout<< len <<std::endl;
    double projNum7[len];
    // double* projNum7 = new double [len];
    for (int i =0; i<len;i++) {
        projNum7[i] = Proj_Value2(p1[i].x,p1[i].y,p1[i].z, c0, step, phiSupport);  
    }

    std::chrono::steady_clock::time_point end6 = std::chrono::steady_clock::now();
    std::cout << "5:Time difference (p0 inside)  projNum_array   = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - begin6 ).count()
    << "[ms]" << std::endl;
    
}


double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport)
{
    
    
    double sum{0};
    int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
    int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
    int zzc = floor(zz), zzf = SampRate * (zz - zzc);
    
     for (int ii = 0; ii < 1e6; ++ii)
    for(int i = 0; i < phiSupport; ++i)
        for(int j = 0; j < phiSupport; ++j)
            for(int k = 0; k < phiSupport; ++k)
            {
                sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                        * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
            }
    return sum;
}

double Proj_Value2(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport)
{
    
    
    double sum{0};
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
    return sum;
}
    

void Proj_Value3( std::vector<Particle>& p0, double* s, std::vector<int> step, int phiSupport,std::vector<double>& result)
{
    
    double sum{0};

    int len = p0.size();

    for(size_t n = 0; n < len; ++n){

    int xxc = floor(p0[n].x), xxf = SampRate * (p0[n].x - xxc);      
    int yyc = floor(p0[n].y), yyf = SampRate * (p0[n].y - yyc);      
    int zzc = floor(p0[n].z), zzf = SampRate * (p0[n].z - zzc);
    
    for(int i = 0; i < phiSupport; ++i)
        for(int j = 0; j < phiSupport; ++j)
            for(int k = 0; k < phiSupport; ++k)
            {
                sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                        * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
            }
    
    result.push_back(sum);
    }

}


void result_interpret1(const double* s, std::vector<Particle>& p0, std::vector<double>& result)
{
    
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const int SampRate   = (phi.size() - 2) / phiSupport;
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates


    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

    
    auto begin5 = std::chrono::steady_clock::now();
    for(Particle n : p0)
    {
        
        double sum {0};

        // rescale the particle coordinates to MRA framework
        double xx = n.x * ScaleFactor;
        double yy = n.y * ScaleFactor;
        double zz = n.z * ScaleFactor;

        // get the 'Coarse' and 'Finer' coodinate, e.g. if
        // position p == 63.25678, then c == 63, f == 256
        int xxc = floor(xx) ,xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy) ,yyf = SampRate * (yy - yyc);      
        int zzc = floor(zz) ,zzf = SampRate * (zz - zzc);

        for(int i = 0; i < phiSupport; ++i)
            
            for(int j = 0; j < phiSupport; ++j)
                
                for(int k = 0; k < phiSupport; ++k)
                {
                    sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                            * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }

        result.push_back(sum);
   
    }
        auto end5 = std::chrono::steady_clock::now();

        std::cout << "Time difference 5 result_intepret(test)  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - begin5).count()
    << "[ms]" << std::endl;
}






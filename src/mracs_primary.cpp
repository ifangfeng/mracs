#include"mracs_primary.hpp"
#include"kernel.hpp"


#define LOWER_RESOLUTION 8

void welcome()
{
    
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "------------------ MRA_CS Starting... -------------" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    
}

void read_parameter()
{
    
    std::string ParamFile {"param.txt"};
    int npri {0};

    std::ifstream iprmfs (ParamFile);
    if(!iprmfs)
    {
        std::cout << "Reading parameter file " + ParamFile + " with error, Abort " << std::endl;
        std::terminate();
    }
    std::string temp, itemp;
    while(iprmfs >> itemp)
    {
        if(itemp=="Resolution")
        {
            iprmfs >> temp;
            Resolution = atoi(temp.c_str());
            npri++;
        }
        else if(itemp=="phiGenus")
        {
            iprmfs >> temp;
            phiGenus = atoi(temp.c_str());
            npri++;
        }
        else if(itemp=="SampRate")
        {
            iprmfs >> temp;
            SampRate = atoi(temp.c_str());
            npri++;
        }
        else if(itemp == "KernelFunc")
        {
            iprmfs >> temp;
            KernelFunc = atoi(temp.c_str());
            npri++;
        }
        else if(itemp == "Radius")
        {
            iprmfs >> temp;
            Radius = atof(temp.c_str());
            npri++;
        }
        else if(itemp == "DIREC")
        {
            iprmfs >> temp;
            DIREC = temp;
            npri++;
        }
        else if(itemp == "MilleCata")
        {
            iprmfs >> temp;
            MilleCata = temp;
            npri++;
        }
        else if(itemp == "SimBoxL")
        {
            iprmfs >> temp;
            SimBoxL  = atof(temp.c_str());
            npri++;
        }
        else if(itemp == "Threads")
        {
            iprmfs >> temp;
            Threads = atoi(temp.c_str());
            npri++;
        }
    }
    
    
    std::cout << "Reading param.txt " << std::endl;
    std::cout << "-> Resolution  =  " << Resolution << std::endl;
    std::cout << "-> phiGenus    =  " << phiGenus << std::endl;
    std::cout << "-> SampRate    =  " << SampRate << std::endl;
    std::cout << "-> KernelFunc  =  " << KernelFunc << std::endl;
    std::cout << "-> Radius      =  " << Radius << std::endl;
    std::cout << "-> DIREC       =  " << DIREC << std::endl;
    std::cout << "-> MilleCata   =  " << MilleCata << std::endl;
    std::cout << "-> SimBoxL     =  " << SimBoxL << std::endl;
    std::cout << "-> Threads     =  " << Threads << std::endl;
    
    omp_set_num_threads(Threads);

    if(npri == PARAMETERS_NUM)     
    {
        std::cout << "Done!"<< std::endl;
    }
    else
    {
        std::cout << "!Missing paramters, Abort "<< std::endl;
        std::terminate();
    }

    const int L {1<<Resolution};

    RESOL = "L" + std::to_string(L);
    RADII = "R" + std::to_string(Radius);
    GENUS = "DaubG" + std::to_string(phiGenus);
    std::cout << "---Grid 3d N = "<< L << "^3"<< std::endl;

}



//=======================================================================================
//|||||||||||||||||||||||||||||| read in wavelets phi data ||||||||||||||||||||||||||||||
//=======================================================================================
//---- numerical value of Daubechies scaling function of genus X, locate in closed
//---- interval [0, 2X-1] with sampling rate 10000 points per unit length (points / 1)
//---------------------------------------------------------------------------------------
std::vector<double> read_in_phi(const int phiGenus, const std::string DIREC)
{
    const int phiStart   {0};                     //Wavelet Phi0() has compact support, 
    const int phiEnd     {2*phiGenus - 1};        //start in x == 0, end in phi_end == 2n-1
    const int phiSupport {phiEnd - phiStart};     //Wavelet Phi0() has compact support 
    
    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;
    
    std::string iwfname {DIREC + "/wavelets/DaubechiesG" + std::to_string(phiGenus) + "Phi.bin"};
    std::ifstream iwfs{iwfname, std::ios_base::binary};
    if(!iwfs) 
    {   
        std::cout << "!Reading " + iwfname + " with error..." << std::endl;
        std::terminate();
    }
    
    double temp;
    std::vector<double> phi;
    while (iwfs.read(as_bytes(temp), sizeof(double)))
    {
        phi.push_back(temp);
    }
    
    phi.push_back(phiStart);
    phi.push_back(phiEnd);

    std::cout << "---Daubechies phi genus: " << phiGenus << std::endl;
    std::cout << "---Wavelet phi supports: [0," << phiSupport << ") \n";
    std::cout << "---Sampling points: " << phi.size() - 2 << std::endl;

    return phi;
}

//=======================================================================================
//||||||||||||||| sfc3d (scaling function coefficients of density field) ||||||||||||||||
//=======================================================================================
double* scaling_function_coefficients(std::vector<Particle>& p, std::vector<double>& phi, const int J , const double SimBoxL)
{
    const int L          = 1 << J;
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const double ScaleFactor {L/SimBoxL};   //used to rescale particle coordinates

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
    {
        step[i] = i * SampRate;
    }

    auto s = new double[L*L*L]();         // density field coefficients in v_j space
    

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    double s_temp[phiSupport * phiSupport * phiSupport];        // openmp seems not allow Array Reductions allocated on heap

    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction (+:s_temp)
    #endif
    
    for(int n = 0; n < p.size(); ++n)
    {
        for(int ii = 0; ii < phiSupport * phiSupport * phiSupport; ++ii) s_temp[ii] = 0;
        // rescale the particle coordinates to MRA framework
        double xx = p[n].x * ScaleFactor;
        double yy = p[n].y * ScaleFactor;
        double zz = p[n].z * ScaleFactor;

        // get the 'Coarse' and 'Finer' coodinate, e.g. if
        // position p == 63.25678, then c == 63, f == 256
        int xxc = xx, xxf = SampRate * (xx - xxc);      
        int yyc = yy, yyf = SampRate * (yy - yyc);      
        int zzc = zz, zzf = SampRate * (zz - zzc);
    
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                {
                    s_temp[i*phiSupport*phiSupport + j*phiSupport + k] += phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                    s[((xxc - i) & (L - 1)) * L * L+ ((yyc - j) & (L - 1)) * L + ((zzc - k) & (L - 1))] += s_temp[i*phiSupport*phiSupport + j*phiSupport + k];        
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 1 sfc3d    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    return s;
}

//read in Real Value Vector v, return an array s, which store the continuous fourier
//transform of v, located in frequency space [k0, k1] with N_k + 1 sampling points
double* Spectrum1(std::vector<double>& v, double k0, double k1, int N_k)
{
    int N_x= v.size()-2;
    double x0 = v[v.size()-2];
    double x1 = v[v.size()-1];
    

    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};

    double Real, Image, Temp;
    double Phase = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    
    double* s = new double[(N_k +1) * 2]();
    for(int i = 0; i <= N_k; ++i)
    {
        for(int j = 0; j < N_x; ++j)
        { 
            Phase = -TWOPI * (k0 + i * Delta_k) * (x0 + j * Delta_x);
            Real  = cos(Phase);
            Image = sin(Phase);
            Temp = v[j] * Delta_x;
            s[2*i]   += Temp * Real;  //- v[2*j+1] * Image;
            s[2*i+1] += Temp * Image; //+ v[2*j+1] * Real;
        }
    }
    return s;
}

//faster but slightly less accuate than Spectrum1()
double* Spectrum(std::vector<double>& v, double k0, double k1, int N_k)
{
    int N_x = v.size()-2;
    double x0 = v[v.size()-2];
    double x1 = v[v.size()-1];
    

    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};
    
    double Real, Image, Temp;
    double Phase = -TWOPI * k0 * x0;
    double PhaseOut = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    const double dPhaseOut = -TWOPI * Delta_k * x0;
    const double dDeltaPhase = -TWOPI * Delta_k * Delta_x;
    

    double* s = new double[(N_k +1) * 2]();
    for(int i = 0; i <= N_k; ++i)
    {
        for(int j = 0; j < N_x; ++j)
        {
            Real  = cos(Phase);
            Image = sin(Phase);
            Temp = v[j] * Delta_x;
            s[2*i]   += Temp * Real;  //- v[2*j+1] * Image;
            s[2*i+1] += Temp * Image; //+ v[2*j+1] * Real;
            Phase += DeltaPhase;
        }
        PhaseOut += dPhaseOut;
        Phase = PhaseOut;
        DeltaPhase += dDeltaPhase;
    }
    return s;
}

//read in Real Value Vector v, calculate its power spectrum located in 
//frequency space [k0, k1] with N_k + 1 sampling points, return as array p[]
double* PowerSpectrum(std::vector<double>& v, double k0, double k1, int N_k)
{
    double* s = Spectrum(v, k0, k1, N_k);
    double* p = new double[N_k + 1];
    for(int i = 0; i <= N_k; ++i)
        p[i] = s[2*i] * s[2*i] + s[2*i+1] * s[2*i+1];
    delete[] s;
    return p;
}

//=======================================================================================
//||||||||||||||| wfc3d (scaling function coefficients of window function) ||||||||||||||
//=======================================================================================
double* window_function_coefficients(std::vector<double>& phi, const int J, const double SimBoxL, const double Radius)
{
    const int L {1<<J};
    const int bandwidth = 1;
    const double DeltaXi = 1./(L);
    const double rescaleR {Radius * L / SimBoxL};

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    double (*WindowFunction)(double, double, double, double){nullptr};

    // WindowFunction point to correct kernel
    if(KernelFunc == 0) 
    {
        WindowFunction = WindowFunction_Shell;
    }
    else if(KernelFunc == 1) 
    {
        WindowFunction = WindowFunction_Sphere;
    }
    else if(KernelFunc == 2) 
    {
        WindowFunction = WindowFunction_Gaussian;
    }


    double* PowerPhi = PowerSpectrum(phi, 0, bandwidth, L * bandwidth);     // Fourier of autocorrelation of WaveletsPhi
    auto WindowArray = new double[(L+1) * (L+1) * (L+1)];                   // temporary array stores kernel's convolution
    auto w           = new double[L * L * L]();                             // window function coefficients in v_j space
    double temp;

    #ifdef IN_PARALLEL
    #pragma omp parallel for private(temp)
    #endif

    for(int i = 0; i <= L; ++i)
        for(int j = 0; j <= L; ++j)
            for(int k = 0; k <= L; ++k)
            {
                temp = 0;
                for(int ii = 0; ii < bandwidth; ++ii)
                    for(int jj = 0; jj < bandwidth; ++jj)
                        for(int kk = 0; kk < bandwidth; ++kk)
                            temp += PowerPhi[ii * L + i] * PowerPhi[jj * L + j] * PowerPhi[kk * L + k]
                                    * WindowFunction(rescaleR, (ii * L + i) * DeltaXi, (jj * L + j) * DeltaXi, (kk * L + k) * DeltaXi);
                
                WindowArray[i * (L+1) * (L+1) + j * (L+1) + k] = temp;
            }

    #ifdef IN_PARALLEL     
    #pragma omp parallel for
    #endif

    for(int i = 0; i < L; ++i)
        for(int j = 0; j < L; ++j)
            for(int k = 0; k < L; ++k)
            {
                w[i * L * L + j * L + k] = WindowArray[i * (L+1) * (L+1) + j * (L+1) + k]
                                           + WindowArray[(L-i) * (L+1) * (L+1) + j * (L+1) + k]
                                           + WindowArray[i * (L+1) * (L+1) + (L-j) * (L+1) + k]
                                           + WindowArray[i * (L+1) * (L+1) + j * (L+1) + (L-k)]
                                           + WindowArray[(L-i) * (L+1) * (L+1) + (L-j) * (L+1) + k]
                                           + WindowArray[(L-i) * (L+1) * (L+1) + j * (L+1) + (L-k)]
                                           + WindowArray[i * (L+1) * (L+1) + (L-j) * (L+1) + (L-k)]
                                           + WindowArray[(L-i) * (L+1) * (L+1) + (L-j) * (L+1) + (L-k)];
            }
    

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()
    << "[ms]" << std::endl;

    delete[] WindowArray;

    return w;
}

//calculate two array's inner product and store result in the first one
void inner_product0(double* v0, double* v1, int N)
{
    for(int i = 0; i < N; ++i)
        v0[i] *= v1[i];
}

//s and w are 3d Real array, s in physical space while w in frequency space
//convol == fftback(inner_product(fft(s), w)), Matrix3D == L*L*L , L == 2^J
//convolution result store in s
void specialized_convolution_3d(double* s, double* w, int J)
{
    const int L {1<<J};
    const int N {L*L*L};

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();

    //auto sc = new double[N*2];
    auto sc = fftw_alloc_complex(N);

    fftw_plan pl1 = fftw_plan_dft_3d(L, L, L, sc, sc, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan pl2 = fftw_plan_dft_3d(L, L, L, sc, sc, FFTW_BACKWARD, FFTW_MEASURE);

    #pragma omp parallel for
    for(int i = 0; i < N; ++i)
        sc[i][0] = s[i];

    
    fftw_execute(pl1);
    //FFT3D_CUBIC(sc, J, 1);
    #pragma omp parallel for
    for(int i = 0; i < N; ++i)
    {
        sc[i][0] *= w[i];
        sc[i][1] *= w[i]; 
    }
    
    fftw_execute(pl2);
    //FFT3D_CUBIC(sc, J, 0);
    #pragma omp parallel for
    for(int i = 0; i < N; ++i)
        s[i] = sc[i][0]/N;

    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3 convl3d  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count()
    << "[ms]" << std::endl;

    fftw_free(sc);
}

//=======================================================================================
//|||||||||||||||||||||||||||||||||||| Result interpret |||||||||||||||||||||||||||||||||
//=======================================================================================
void result_interpret(const double* s, const int J, const double SimBoxL, std::vector<double>& phi,
                     std::vector<Particle>& p0, std::vector<double>& result, const double total)
{
    const int L          = 1 << J;
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const int SampRate   = (phi.size() - 2) / phiSupport;
    const double ScaleFactor {L/SimBoxL};   //used to rescale particle coordinates


    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

    const int N {L*L*L};

    
    for(int n = 0; n < p0.size(); ++n)
    {
        double sum {0};

        // rescale the particle coordinates to MRA framework
        double xx = p0[n].x * ScaleFactor;
        double yy = p0[n].y * ScaleFactor;
        double zz = p0[n].z * ScaleFactor;

        // get the 'Coarse' and 'Finer' coodinate, e.g. if
        // position p == 63.25678, then c == 63, f == 256
        int xxc = xx, xxf = SampRate * (xx - xxc);      
        int yyc = yy, yyf = SampRate * (yy - yyc);      
        int zzc = zz, zzf = SampRate * (zz - zzc);

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                {
                    sum += s[((xxc - i) & (L - 1)) * L * L+ ((yyc - j) & (L - 1)) * L + ((zzc - k) & (L - 1))]
                            * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }

        result.push_back(sum);
        
    }

    //std::cout << "Check: " << std::endl;

    //double Expect, Calcul{0};
    //Expect = 4./3. * M_PI * pow((Radius * ScaleFactor), 3) / N * total;
    //std::cout << "Expectation of convolved Galaxy Numbers: " << Expect << std::endl;

    //for(auto x : result) Calcul += x;
    //Calcul /= result.size();
    //std::cout << "Calculation of convolved Galaxy Numbers: " << Calcul << std::endl;


}


//calculate index of box that must be inside the sphere R, and that might be intersect with sphere boundray
void fill_index_set(const double R, std::vector<Index>& inner_index, std::vector<Index>& cross_index)
{
    const int M = R + 1;
    for(int i = -M; i <= M; ++i)
        for(int j = -M; j <= M; ++j)
            for(int k = -M; k <= M; ++k)
            {
                if(sqrt(i * i + j * j + k * k) <= R - sqrt(3.))
                {
                    inner_index.push_back(Index{i, j, k});
                }
                else if(sqrt(i * i + j * j + k * k) < R + sqrt(3.))
                {
                    cross_index.push_back(Index{i, j, k});
                }
            }
}


//particle periodic in box size L, calculate the number of particles in sphere center at p0 with radius R
void count_in_sphere(const double R, const double SimBoxL, std::vector<Particle>& p, std::vector<Particle>& p0, std::vector<double>& result)
{
    std::vector<Index> inner_index;
    std::vector<Index> cross_index;
    std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();

    fill_index_set(R/SimBoxL, inner_index, cross_index);

    double count[p0.size()];
    double temp;

    #ifdef IN_PARALLEL
    #pragma omp parallel for private(temp)
    #endif

    for(size_t n = 0; n < p0.size(); ++n)
    {
        temp = 0;
        for(size_t i = 0; i < p.size(); ++i)
            for(size_t m = 0; m < cross_index.size(); ++m)
            {
                double xx = p[i].x + cross_index[m].i * SimBoxL - p0[n].x;
                double yy = p[i].y + cross_index[m].j * SimBoxL - p0[n].y;
                double zz = p[i].z + cross_index[m].k * SimBoxL - p0[n].z;
                if((abs(xx) < R) && (abs(yy) < R) && (abs(zz) < R))
                    if(xx*xx + yy*yy + zz*zz < R*R)
                        ++temp;
            }
        count[n]= temp + inner_index.size() * p.size();
    }
    for(int i = 0; i < p0.size(); ++i)
        result.push_back(count[i]);

    std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 count    = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;
}


void stupid_count(const double R, std::vector<Particle>& p, std::vector<Particle>& p0, std::vector<double>& stupid_result)
{
    double count[p0.size()];

    #ifdef IN_PARALLEL
    #pragma omp parallel for 
    #endif

    for(size_t n = 0; n < p0.size(); ++n)
    {
        count[n] = 0;
        for(size_t i = 0; i < p.size(); ++i)
        {
            double xx = p[i].x - p0[n].x;
            double yy = p[i].y - p0[n].y;
            double zz = p[i].z - p0[n].z;
            if(xx*xx + yy*yy + zz*zz < R*R)
                ++count[n];
        }
    }

    for(size_t i = 0; i < p0.size(); ++i)
        stupid_result.push_back(count[i]);

}














/*
void output_result(const double* v, const int J, )
{
    
    std::string ofname = DIREC + "/output/sfc/DensityField" + GENUS + RESOL + ".bin";;
    std::ofstream ofs(ofname);
    if(!ofs)
    {
        std::cout << "!Writting sfc file with error, Abort " << std::endl;
        std::terminate();
    }
    
}

    std::cout << "============direct count test============" << '\n';
    std::cout << "rescale_R := R/SimBoxL= " << R/SimBoxL << '\n';
    std::cout << "------------inner_index------------------" << '\n';
    std::cout << "TotalNum= " << inner_index.size() << '\n';

    for(auto x : inner_index) std::cout << x.i << " " << x.j << " " << x.k << '\n';

    std::cout << "------------cross_index------------------" << '\n';
    std::cout << "TotalNum= " << cross_index.size() << '\n';
    for(auto x : cross_index) std::cout << x.i << " " << x.j << " " << x.k << '\n';
*/

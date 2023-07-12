#include"MRACS_Main.h"
#include"MRACS_Kernel.h"
#include"MRACS_Readin.h"



double* PowerPhi=nullptr;
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
        else if(itemp=="BaseType")
        {
            iprmfs >> temp;
            BaseType = atoi(temp.c_str());
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
        else if(itemp == "DataDirec")
        {
            iprmfs >> temp;
            DataDirec = temp;
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
    std::string BaseType_String = !BaseType ? "B_Spline" : "Daubechies";
    
    std::cout << "Reading param.txt " << std::endl;
    std::cout << "-> Resolution  =  " << Resolution << std::endl;
    std::cout << "-> BaseType    =  " << BaseType_String << std::endl;
    std::cout << "-> phiGenus    =  " << phiGenus << std::endl;
    std::cout << "-> SampRate    =  " << SampRate << std::endl;
    std::cout << "-> KernelFunc  =  " << KernelFunc << std::endl;
    std::cout << "-> Radius      =  " << Radius << std::endl;
    std::cout << "-> DIREC       =  " << DIREC << std::endl;
    std::cout << "-> DataDirec   =  " << DataDirec << std::endl;
    std::cout << "-> SimBoxL     =  " << SimBoxL << std::endl;
    std::cout << "-> Threads     =  " << Threads << std::endl;
    
    omp_set_num_threads(Threads);

    if(npri == PARAMETERS_NUM)     
    {
        std::cout << "Initializing..."<< std::endl;
    }
    else
    {
        std::cout << "!Missing paramters, Abort "<< std::endl;
        std::terminate();
    }


    GridLen = 1 << Resolution;
    GridVol = 1UL << Resolution*3;

    RESOL = "L" + std::to_string(GridLen);
    RADII = "R" + std::to_string(Radius);
    GENUS = "DaubG" + std::to_string(phiGenus);
    
    if(BaseType == 0)
    {
        phi = B_Spline(phiGenus, SampRate);
        phiSupport = phiGenus + 1;
    }
    else if (BaseType == 1)
    {
        phi = Daubechies_Phi(phiGenus);
        phiSupport = 2*phiGenus - 1;
    }

    PowerPhi = PowerPhiFunc(GridLen);

    std::cout << "Done!"<< std::endl;
}



//=======================================================================================
//|||||||||||||||||||||||||||||| B-Spline phi data ||||||||||||||||||||||||||||||
//=======================================================================================
//---- numerical value of B-Spline scaling function of genus n, locate in closed interval
//---- [-(n+1)/2, (n+1)/2] with sampling rate 1000 points per unit length (points / 1)
//---------------------------------------------------------------------------------------
std::vector<double> B_Spline(const int n, const int sampleRate)
{
    double* spline = new double[(n+1)*sampleRate+1];
    const double delta_x = 1./sampleRate;

    for(int i = 1; i < sampleRate; ++i)
    {
        spline[i] = 1.;
    }
    spline[0] = 0.5;
    spline[sampleRate] = 0.5;

    for(int m = 1; m <= n; ++m)
    {
        for(size_t i = 0; i <= m * sampleRate; ++i)
        {
            spline[i] *= i * delta_x / m;
        }
        for(size_t i = sampleRate; i <= (m+1) * sampleRate / 2; ++i)
        {
            spline[i] += spline[(m+1) * sampleRate - i];
        }
        for(size_t i = 0; i <= (m+1) * sampleRate / 2; ++i)
        {
            spline[(m+1) * sampleRate - i] = spline[i];
        }
    }
    std::vector<double> phi;
    for(size_t i = 0; i <= (n+1) * sampleRate; ++i)
    {
        phi.push_back(spline[i]);
    }
    delete[] spline;

    return phi;
}


//=======================================================================================
//||||||||||||||||||||||||||||| Daubechies phi data |||||||||||||||||||||||||||||
//=======================================================================================
//---- numerical value of Daubechies scaling function of genus n, locate in closed
//---- interval [0, 2n-1] with sampling rate 1000 points per unit length (points / 1)
//---------------------------------------------------------------------------------------
std::vector<double> Daubechies_Phi(const int phiGenus)
{
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
    
    return phi;
}

//=======================================================================================
// sfc of unweighted grid point
//=======================================================================================
double* sfc_grid_coordinate(std::vector<int64_t>& ps)
{   
    auto s = new double[GridVol]();         // density field coefficients in v_j space
    
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    double s_temp[phiSupport * phiSupport * phiSupport];
    

    #pragma omp parallel for reduction (+:s_temp)
    for(size_t n = 0; n < ps.size(); ++n)
    {
        for(int ii = 0; ii < phiSupport * phiSupport * phiSupport; ++ii) s_temp[ii] = 0;
        int64_t xx,yy,zz;
        zz = ps[n]&(GridLen-1);
        yy = (ps[n]&((GridLen-1)<<Resolution))>>Resolution;
        xx = (ps[n]&((GridLen-1)<<(Resolution*2)))>>(Resolution*2);

        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                    s_temp[i*phiSupport*phiSupport + j*phiSupport + k] += phi[i * SampRate] * phi[j * SampRate] * phi[k * SampRate];
    
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                    s[((xx - i)&(GridLen - 1))*GridLen*GridLen + ((yy - j)&(GridLen - 1))*GridLen + ((zz - k)&(GridLen - 1))] +=
                    s_temp[i*phiSupport*phiSupport + j*phiSupport + k];        
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference grid_sfc   = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    return s;
}

//=======================================================================================
//||||||||||||||| sfc3d (scaling function coefficients of density field) ||||||||||||||||
//=======================================================================================
double* sfc_offset(std::vector<Particle>& p, Offset v)
{   
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates
    
    // double total {0};
    // #pragma omp parallel for reduction (+:total)
    // for(auto x : p) total += x.weight;
    // total /= p.size();

    auto s = new double[GridVol]();         // density field coefficients in v_j space

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    double s_temp[phiSupport * phiSupport * phiSupport];

    #pragma omp parallel for reduction (+:s_temp)
    for(size_t n = 0; n < p.size(); ++n)
    {
        for(int ii = 0; ii < phiSupport * phiSupport * phiSupport; ++ii) s_temp[ii] = 0;
        // rescale the particle coordinates to MRA framework
        double xx = (p[n].x + v.dx) * ScaleFactor;
        double yy = (p[n].y + v.dy) * ScaleFactor;
        double zz = (p[n].z + v.dz) * ScaleFactor;

        // get the 'Coarse' and 'Finer' coodinate, e.g. if
        // position p == 63.25678, then c == 63, f == 256
        int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
        int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
        int zzc = floor(zz), zzf = SampRate * (zz - zzc);
    
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                {
                    s_temp[i*phiSupport*phiSupport + j*phiSupport + k] += 
                    phi[xxf + i * SampRate] * phi[yyf + j * SampRate] * phi[zzf + k * SampRate] * p[n].weight;// / total;
                }
        for(int i = 0; i < phiSupport; ++i)
            for(int j = 0; j < phiSupport; ++j)
                for(int k = 0; k < phiSupport; ++k)
                    s[((xxc - i)&(GridLen - 1))*GridLen*GridLen + ((yyc - j)&(GridLen - 1))*GridLen + ((zzc - k)&(GridLen - 1))] +=
                    s_temp[i*phiSupport*phiSupport + j*phiSupport + k];        
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 1 sfc3d    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    return s;
}

double* sfc(std::vector<Particle>& p)
{
    Offset a(0., 0., 0.);
    return sfc_offset(p, a);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// read in Real-Value vector v, return an array s, which store the continuous fourier
// transform of v, located in frequency space [k0, k1] with N_k + 1 sampling points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double* Spectrum1(std::vector<double>& v, double k0, double k1, size_t N_k)
{
    size_t N_x = v.size();
    double x0 = 0;
    double x1 = phiSupport;
    

    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};

    double Real, Image, Temp;
    double Phase = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    
    double* s = new double[(N_k +1) * 2]();
    for(size_t i = 0; i <= N_k; ++i)
    {
        for(size_t j = 0; j < N_x; ++j)
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// faster but slightly less accuate than Spectrum1()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double* Spectrum(std::vector<double>& v, double k0, double k1, size_t N_k)
{
    size_t N_x = v.size();
    double x0 = 0;
    double x1 = phiSupport;
    

    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};
    
    double Real, Image, Temp;
    double Phase = -TWOPI * k0 * x0;
    double PhaseOut = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    const double dPhaseOut = -TWOPI * Delta_k * x0;
    const double dDeltaPhase = -TWOPI * Delta_k * Delta_x;
    

    double* s = new double[(N_k +1) * 2]();
    for(size_t i = 0; i <= N_k; ++i)
    {
        for(size_t j = 0; j < N_x; ++j)
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

//=======================================================================================
// read in Real Value Vector v, calculate its power spectrum located in 
// frequency space [k0, k1] with N_k + 1 sampling points, return as array p[]
//=======================================================================================
double* PowerSpectrum(std::vector<double>& v, double k0, double k1, size_t N_k)
{
    double* s = Spectrum(v, k0, k1, N_k);
    double* p = new double[N_k + 1];
    for(size_t i = 0; i <= N_k; ++i)
        p[i] = s[2*i] * s[2*i] + s[2*i+1] * s[2*i+1];
    delete[] s;
    return p;
}




//=======================================================================================
// window array (convolution kernel)
//=======================================================================================
double* windowArray(const double Radius, const double theta)
{
    const double DeltaXi = 1./GridLen;
    const double RGrid {Radius * GridLen/SimBoxL};

    auto WindowArray = new double[(GridLen+1) * (GridLen+1) * (GridLen+1)];

    if(KernelFunc <= 2)
    {
        double (*WindowFunction)(double, double, double, double){nullptr};
        if(KernelFunc == 0) WindowFunction = WindowFunction_Shell;
        else if(KernelFunc == 1) WindowFunction = WindowFunction_Sphere;
        else if(KernelFunc == 2) WindowFunction = WindowFunction_Gaussian;
    
        #pragma omp parallel for
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = WindowFunction(RGrid, i * DeltaXi, j * DeltaXi, k * DeltaXi);
    }
    else if(KernelFunc == 3)
    {
        double fz[GridLen+1];
        double dXitwo{pow(DeltaXi,2)};
        for(size_t i = 0; i <= GridLen; ++i) fz[i] = cos(TWOPI * RGrid * cos(theta) * i*DeltaXi);
        double* fxy = new double[(GridLen+1) * (GridLen+1)];

        #pragma omp parallel for 
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                fxy[i * (GridLen+1) + j] = std::cyl_bessel_j(0,TWOPI*sin(theta)*RGrid*sqrt(i*i*dXitwo+j*j*dXitwo));

        #pragma omp parallel for 
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = fxy[i * (GridLen+1) + j] * fz[k]; 
        delete[] fxy;
    }
    else if(KernelFunc == 4)
    {
        double fz[GridLen+1];
        double dXitwo{pow(DeltaXi,2)};
        const double HGrid {theta * GridLen/SimBoxL};
        for(size_t i = 0; i <= GridLen; ++i) fz[i] = sin(M_PI*i*DeltaXi*HGrid)/(M_PI*i*DeltaXi*HGrid);fz[0]=1;
        double* fxy = new double[(GridLen+1) * (GridLen+1)];

        #pragma omp parallel for 
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                fxy[i * (GridLen+1) + j] = std::cyl_bessel_j(1,TWOPI*RGrid*sqrt(i*i*dXitwo+j*j*dXitwo))/(M_PI*RGrid*sqrt(i*i*dXitwo+j*j*dXitwo));fxy[0]=1;
        
        #pragma omp parallel for 
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = fxy[i * (GridLen+1) + j] * fz[k];
        delete[] fxy;
    }

    return WindowArray;
}

// Power of window function W(Radius,theta)
double *window_Pk(const double Radius, const double theta)
{
    auto winPk = windowArray(Radius,theta);
    #pragma omp paraller for
    for(size_t i = 0; i < (GridLen+1)*(GridLen+1)*(GridLen+1); ++i) winPk[i] = pow(winPk[i], 2);
    
    return winPk;
}

//=======================================================================================
// calculate m_th B_Spline dual's power spectrum located in frequency
// space [k0, k1] with N_k + 1 sampling points, return as array p[]
//=======================================================================================
double* B_Spline_Dual_Power_Spectrum(double m, double k0, double k1, size_t N_k)
{
    double* p = new double[N_k + 1];
    double* a = new double[N_k + 1];
    
    std::vector<double> c = B_Spline(2*m+1, 1);
    for(size_t i = 0; i <= N_k; ++i)
    {
        a[i] = c[m+1];
    }
    const double Delta_k {(k1-k0)/N_k};
    double k = k0;
    for(int n = 1; n <= m; ++n)
    {
        for(size_t i = 0; i <= N_k; ++i)
        {
            a[i] += 2*c[m+n+1] * cos(TWOPI*n*k);
            k += Delta_k;
        }
        k = k0;
    }
    k = k0;
    for(size_t i = 0; i <= N_k; ++i)
    {
        if(k == 0)
        {
            p[i] = pow(a[i], -2);
        }
        else
        {
            p[i] = pow(sin(k*M_PI)/(k*M_PI), 2*(m+1))/pow(a[i], 2);
        }
        k += Delta_k;
    }
    
    delete[] a;
    return p;
}

// power_phi array
double* PowerPhiFunc(const size_t N)
{
    double* PowerSpectra = nullptr;
    if(BaseType == 0)
    {
        PowerSpectra = B_Spline_Dual_Power_Spectrum(phiGenus, 0, 1, N);
    }
    else if (BaseType == 1)
    {
        PowerSpectra = PowerSpectrum(phi, 0, 1, N);
    }

    auto a = new double[(N+1)*(N+1)*(N+1)];

    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i <= N; ++i)
        for(size_t j = 0; j <= N; ++j)
            for(size_t k = 0; k <= N; ++k)
            {
                a[i * (N+1) * (N+1) + j * (N+1) + k] = PowerSpectra[i] * PowerSpectra[j] * PowerSpectra[k];
            }
    return a;
}


void force_resoluton_J(int j)
{
    if(Resolution != j){
        Resolution = j;
        GridLen = 1 << Resolution;
        GridVol = 1UL << Resolution*3;
        delete[] PowerPhi;
        PowerPhi = PowerPhiFunc(GridLen);
        std::cout << "!MRACS resolution has been forced to " << j << "\n";
    }
}

// 0-Shell; 1-Sphere; 2-Gaussian; 3-DualRing
void force_kernel_type(int x)
{
    if(x != KernelFunc)
    {
        KernelFunc = x;
        delete[] PowerPhi;
        PowerPhi = PowerPhiFunc(GridLen);
        std::cout << "!kernel function has been forced to " << x << "\n";
    }

}

// 'a' specify BaseType{0:BSpline,1:DaubWavelet} while 'n' phiGenus
void force_base_type(int a, int n)
{
    if(a < 0 || a > 1 || n < 0)
    {
        std::cout << "!Illegal input, a should be '0' or '1' and n is a non-negative integer" << std::endl;
    }
    else if(a != BaseType || n != phiGenus )
    {
        BaseType = a;
        phiGenus = n;
        std::string BaseType_String = !BaseType ? "B_Spline" : "Daubechies";
        if(BaseType == 0)
        {
            phi = B_Spline(phiGenus, SampRate);
            phiSupport = phiGenus + 1;
        }
        else if (BaseType == 1)
        {
            phi = Daubechies_Phi(phiGenus);
        }
        delete[] PowerPhi;
        PowerPhi = PowerPhiFunc(GridLen);
        std::cout << "!BaseType has been forced to: " << BaseType_String << n << "\n";
    }
}

double* symmetryFold_lean(double* wA)
{
    auto w = new double[GridLen * GridLen * (GridLen/2+1)]();

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                w[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k]
                = wA[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k]
                + wA[(GridLen-i) * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k]
                + wA[i * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + k]
                + wA[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + (GridLen-k)]
                + wA[(GridLen-i) * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + k]
                + wA[(GridLen-i) * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + (GridLen-k)]
                + wA[i * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + (GridLen-k)]
                + wA[(GridLen-i) * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + (GridLen-k)];
            }

    delete[] wA;
    return w;
}

//=======================================================================================
//||||||||||||||| Fourier transform of window function (without PowerPhi)||||||||||||||
// when cylinder kernel is adopted, 'theta' should be interpret as shape parameter 'H' 
//=======================================================================================
double* wft(const double Radius, const double theta)
{
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto WindowArray = windowArray(Radius, theta);                          
    auto w = symmetryFold_lean(WindowArray);
    
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wft3d    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    return w;
}

//=======================================================================================
//||||||||||||||| wfc3d (scaling function coefficients of window function) ||||||||||||||
// when cylinder kernel is adopted, 'theta' should be interpret as shape parameter 'H' 
//=======================================================================================
double* wfc(const double Radius, const double theta)
{
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto WindowArray = windowArray(Radius, theta);                          

    #pragma omp parallel for 
    for(size_t i = 0; i < (GridLen+1)*(GridLen+1)*(GridLen+1); ++i) WindowArray[i] *= PowerPhi[i];
  
    auto w = symmetryFold_lean(WindowArray);

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d    = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    return w;
}

// calculate two array's inner product return as double
double inner_product(double* v0, double* v1, size_t N)
{
    double sum {0};

    #pragma omp parallel for reduction (+:sum)
    for(size_t i = 0; i < N; ++i) sum += v0[i] * v1[i];

    return sum;
}

// retrun sum of all elements of array w, w has length N
double array_sum(double* w, size_t N)
{
    double sum {0};

    #pragma omp parallel for reduction (+:sum)
    for(size_t i = 0; i < N; ++i) sum += w[i];

    return sum;
}

fftw_complex* sfc_r2c(double* s, bool DELETE_S)
{
    if(!fftw_init_threads())
    {
        std::cout << "fftw_init_threads() with error, Abort" << std::endl;
        std::terminate();
    }

    auto sc = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));

    fftwf_plan_with_nthreads(omp_get_max_threads());

    fftw_plan pl = fftw_plan_dft_r2c_3d(GridLen, GridLen, GridLen, s, sc, FFTW_MEASURE);
    fftw_execute(pl);

    if(DELETE_S) delete[] s;

    return sc;
}

double* convol_c2r(fftw_complex* sc, double* w)
{
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto c = new double[GridVol];

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    {
        sc1[i][0] = w[i] * sc[i][0];
        sc1[i][1] = w[i] * sc[i][1]; 
    }
    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan pl = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, c, FFTW_MEASURE);

    fftw_execute(pl);
    #pragma omp parallel for
    for(size_t i = 0; i < GridVol; ++i) c[i] /= GridVol;

    fftw_free(sc1);
    return c;
}

//=======================================================================================
// s and w are 3d Real array, s in physical space while w in frequency space, return 
// as double* c , convol==fftback(inner_product(fft(s), w)), Matrix3D==L*L*L , L==2^J
//=======================================================================================
double* convol3d(double* s, double* w, bool DELETE_S)
{
    auto begin3 = std::chrono::steady_clock::now();

    auto sc = sfc_r2c(s,false);
    auto c = convol_c2r(sc, w);

    auto end3 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3 convl3d  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count()
    << "[ms]" << std::endl;

    fftw_free(sc);
    if(DELETE_S) delete[] s;

    return c;
}


// given scaling function coefficients s, return the value it represented at position p0
double* project_value(const double* s, std::vector<Particle>& p0, bool DELETE_S)
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
                            * phi[xxf + i * SampRate] * phi[yyf + j * SampRate] * phi[zzf + k * SampRate];
                }

        result[n] = sum;
        sum = 0;
    }

    auto end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 interpret = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    if(DELETE_S) delete[] s;

    return result;
}


// projected value at each grid point
double* prj_grid(const double* s, bool DELETE_S)
{
    auto a = new double[GridVol];

    double sum = 0;
    #pragma omp parallel for reduction (+:sum)
    for(int x = 0; x < GridLen; ++x)
        for(int y = 0; y < GridLen; ++y)
            for(int z = 0; z < GridLen; ++z)
            {
                for(int i = 0; i < phiSupport; ++i)
                    for(int j = 0; j < phiSupport; ++j)
                        for(int k = 0; k < phiSupport; ++k)
                        {
                            sum += s[((x-i) & (GridLen-1)) * GridLen * GridLen + ((y-j) & (GridLen-1)) * GridLen + ((z-k) & (GridLen-1))]
                                    *phi[i * SampRate] * phi[j * SampRate] * phi[k * SampRate];
                        }
                a[x * GridLen * GridLen + y * GridLen + z] = sum;
                sum = 0;
            }
            
    if(DELETE_S) delete[] s;

    return a;
}


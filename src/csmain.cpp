#include"csmain.h"



#define LOWER_RESOLUTION 8

double* PowerPhi=nullptr;
void welcome()
{
    
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "------------------ MRA_CS Starting... -------------" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    
}
std::vector<double> log_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT)
{
    if(Rmin <= 0) {std::cout << "input error, Rmin should be positive number\n";std::terminate();}
    std::vector<double> r(Npt);
    for(int i = 0; i < Npt; ++i) r[i] = Rmin * pow(Rmax/Rmin, static_cast<double>(i)/Npt);
    if(ENDPOINT) r.push_back(Rmax);

    return r;
}

std::vector<double> linear_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT)
{
    std::vector<double> r(Npt);
    for(int i = 0; i < Npt; ++i) r[i] = Rmin + i * (Rmax - Rmin) / Npt;
    if(ENDPOINT) r.push_back(Rmax);

    return r;
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
    
    double total {0};

    #pragma omp parallel for reduction (+:total)
    for(auto x : p) total += x.weight;
    total /= p.size();

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
                    phi[xxf + i * SampRate] * phi[yyf + j * SampRate] * phi[zzf + k * SampRate] * p[n].weight / total;
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
// window array
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

// return the covariance or variance of density field smoothed by a window w(Radius,theta),
// pk_plus is returned by densityCovarianceArray or densityVarianceArray, WinPk is returned by window_Pk
// when covarianceArray is the case, covar_CombinewithKernel() returns the covariance.
double covar_CombinewithKernel(double* pk_plus, double* WinPk)
{
    double sum{0};
    #pragma omp parallel for reduction (+: sum)
    for(int64_t i = 1 - GridLen; i < GridLen; ++i)
        for(int64_t j = 1 - GridLen; j < GridLen; ++j)
            for(int64_t k = 1 - GridLen; k < GridLen; ++k)
            {
                sum += pk_plus[(i+GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (j+GridLen-1) * (2*GridLen-1) + (k+GridLen-1)] *
                WinPk[abs(i) * (GridLen+1) * (GridLen+1) + abs(j) * (GridLen+1) + abs(k)];
            }
    double temp = pk_plus[(GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (GridLen-1) * (2*GridLen-1) + (GridLen-1)] * WinPk[0];

    return sum/temp-1;
}

// Power of window function W(Radius,theta)
double *window_Pk(const double Radius, const double theta)
{
    auto winPk = windowArray(Radius,theta);
    #pragma omp paraller for
    for(size_t i = 0; i < (GridLen+1)*(GridLen+1)*(GridLen+1); ++i) winPk[i] = pow(winPk[i], 2);
    
    return winPk;
}

// FINE is integral finer parameter and takes value from set {0,1,2,3,4,5}
// notice the Maximun finer parameter is limited to 5 
double Pk_variance_2dRH(double* Pk, const double Radius, const double theta, int FINE)
{
    if(FINE < 0 || 5 < FINE){
        FINE = 0;
        std::cout << "illegal finer parameter, reset as default FINE=0 " << std::endl;
    }
    const int STEP = 1 << FINE;
    const int64_t LEN {(GridLen/2+1)*STEP};
    double fz[LEN], fr[LEN];
    double* Pk_map = new double[LEN*LEN];

    if(KernelFunc == 3){
        for(size_t i = 0; i < LEN; ++i) fz[i] = cos(TWOPI * i * Radius/SimBoxL * cos(theta)/STEP);
        for(size_t i = 0; i < LEN; ++i) fr[i] = std::cyl_bessel_j(0,TWOPI*i*sin(theta)*Radius/SimBoxL/STEP);
    }
    else if(KernelFunc == 4){
        for(size_t i = 0; i < LEN; ++i) fz[i] = sin(M_PI*i*theta/SimBoxL)/(M_PI*i*theta/SimBoxL);fz[0]=1;
        for(size_t i = 0; i < LEN; ++i) fr[i] = std::cyl_bessel_j(1,TWOPI*Radius/SimBoxL*i)/(M_PI*Radius/SimBoxL*i);fr[0]=1;
    }
    else{
        std::cout << "Pk_variance_2dRH: KernelFunc == 3 or. 4 only" << std::endl;
        return 0;
    }
    for(size_t i = 0; i < LEN; ++i)
        for(size_t j = 0; j < LEN; ++j){
            double k = sqrt(i*i + j*j);
            int kc = k;
            double kf = k-kc;
            int index = k/STEP;
            double xl = Pk[index];
            double xr = Pk[index +1];
            double pk = (xr - xl)*kf + xl;
            Pk_map[i*LEN + j] = pk;
        }
    double var{0};
    for(int i = 0; i < LEN; ++i)
        for(int j = 1-LEN; j < LEN; ++j){
            double ii = i;
            var += fr[i]*fz[abs(j)]*Pk_map[i*LEN+abs(j)]*TWOPI*ii/pow(STEP,3);
        }

    return var;    
}



// given the sfc array, caluculate its raw covariance array i.e. no window smoothing
double* densityCovarianceArray(fftw_complex* sc1,fftw_complex* sc2)
{
    auto Pk_array = new double[GridVol];
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]
                + sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
            }
    
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k)
            Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
            Pk_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k];

    auto Pk_plus = new double[(2*GridLen-1)*(2*GridLen-1)*(2*GridLen-1)];

    #pragma omp parallel for
    for(int64_t i = 1 - GridLen; i < GridLen; ++i)
        for(int64_t j = 1 - GridLen; j < GridLen; ++j)
            for(int64_t k = 1 - GridLen; k < GridLen; ++k)
            {
                Pk_plus[(i+GridLen-1) * (2*GridLen-1) * (2*GridLen-1) + (j+GridLen-1) * (2*GridLen-1) + (k+GridLen-1)] = 
                Pk_array[((GridLen+i)%GridLen) * GridLen * GridLen + ((GridLen+j)%GridLen) * GridLen + ((GridLen+k)%GridLen)] * 
                PowerPhi[abs(i) * (GridLen+1) * (GridLen+1) + abs(j) * (GridLen+1) + abs(k)];
            }
    
    return Pk_plus;
}

// given the sfc array, caluculate its raw variance array i.e. no window smoothing
double* densityVarianceArray(fftw_complex* sc)
{
    return densityCovarianceArray(sc,sc);
}


//=======================================================================================
// particles first assigned to grid using different window function, then we can take
// advantage of FFT to get its Fourier Coefficiencs and average all orientation
// to have P(k) as function of scalar module k, notice that sfc function do the 
// assignment step exactly
//=======================================================================================
double* densityPowerFFT(fftw_complex* sc)
{
    const double npart = sc[0][0];
    const int64_t PkASize {(GridLen/2 + 1) * GridLen * GridLen};
    double* Pk_array = new double[PkASize];

    #pragma omp parallel for
    for(size_t i = 0; i < PkASize; ++i) Pk_array[i] = pow(sc[i][0], 2) + pow(sc[i][1], 2);

    int64_t klen = GridLen;
    int nk[klen]{0};
    double* Pk = new double[klen]();

    for(size_t i = 0; i < GridLen ; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                int kM = sqrt(i * i + j * j + k * k);
                if(kM < klen){
                Pk[kM] += Pk_array[i * GridLen * (GridLen/2+1) + j * (GridLen/2+1) + k];
                nk[kM] += 1;}
            }
    delete[] Pk_array;

    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(npart,2);
        //Pk[i] -= 1./pow(npart,1); // poission shot noise
    }
    Pk[0]=0;

    return Pk;
}

//=======================================================================================
// calculate density power spectrum using projected density fileds mathematically
//=======================================================================================
double* densityPowerDWT(fftw_complex* sc)
{
    const double npart = sc[0][0];
    double* Pk_array = new double[(GridLen/2 + 1) * GridLen * GridLen];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                Pk_array[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k] = PowerPhi[i * (GridLen + 1) * (GridLen + 1) + j * (GridLen + 1) + k] *
                (pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2));
            }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Pk power, T = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    int64_t klen = GridLen;
    int nk[klen]{0};
    double* Pk = new double[klen]();
    
    for(size_t i = 0; i < GridLen ; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                int kM = sqrt(i * i + j * j + k * k);
                if(kM < klen)
                {
                Pk[kM] += Pk_array[i * GridLen * (GridLen/2+1) + j * (GridLen/2+1) + k];
                nk[kM] += 1;}
            }
    
    delete[] Pk_array;

    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(npart,2);
        
    }
    Pk[0]=0;
   
    return Pk;
}
//Pk[i] -= 1./pow(npart,1); // poission shot noise

//=======================================================================================
// calculate cross power spectrum of sc1 and sc2
//=======================================================================================
double* crossPowerDWT(fftw_complex* sc1, fftw_complex* sc2)
{
    const double npartsq = sc1[0][0] * sc2[0][0];
    double* Pk_array = new double[(GridLen/2 + 1) * GridLen * GridLen];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                Pk_array[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k] = PowerPhi[i * (GridLen + 1) * (GridLen + 1) + j * (GridLen + 1) + k] *
                sqrt((pow(sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + pow(sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2)) *
                (pow(sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + pow(sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2)));
            }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Pk power, T = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    int64_t klen = GridLen;
    int nk[klen]{0};
    double* Pk = new double[klen]();
    
    for(size_t i = 0; i < GridLen ; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                int kM = sqrt(i * i + j * j + k * k);
                if(kM < klen)
                {
                Pk[kM] += Pk_array[i * GridLen * (GridLen/2+1) + j * (GridLen/2+1) + k];
                nk[kM] += 1;}
            }
    
    delete[] Pk_array;

    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= npartsq;
        
    }
    Pk[0]=0;
   
    return Pk;
}

//=======================================================================================
// calculate cross-correlation function c(k) of two density fileds in fourier space,
// using mass assignment and FFT method. where k is a scalar.
//=======================================================================================
double* densityCorrelationFFT(fftw_complex* sc1, fftw_complex* sc2)
{
    auto cross_array = new double[GridVol][2];

    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                cross_array[i * GridLen * GridLen + j * GridLen + k][0] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]+ 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
                cross_array[i * GridLen * GridLen + j * GridLen + k][1] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]-
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
            }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k){
                cross_array[i * GridLen * GridLen + j * GridLen + k][0] = 
                cross_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k][0];
                cross_array[i * GridLen * GridLen + j * GridLen + k][1] = 
                cross_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k][1];
            }
    int klen = GridLen*sqrt(3.);
    int nk[klen];
    auto ccf = new double[klen];
    for(int i = 0; i < klen; ++i)
    {   
        nk[i] = 0;
        ccf[i] = 0;
    }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k){
                int kM = sqrt(i * i + j * j + k * k);
                ccf[kM] += sqrt(pow(cross_array[i * GridLen * GridLen + j * GridLen + k][0],2)+ 
                                pow(cross_array[i * GridLen * GridLen + j * GridLen + k][1],2));
                nk[kM] += 1;
            }
    delete[] cross_array;
    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        ccf[i] /= nk[i];
        ccf[i] /= pow(GridVol,2);
    }
    
    return ccf;
}


//=======================================================================================
// calculate cross-correlation function c(k) of two density fileds in fourier space,
// using projected density fileds mathematically, where k is a scalar.
//=======================================================================================
double* densityCorrelationDWT(fftw_complex* sc1, fftw_complex* sc2)
{
    auto cross_array = new double[GridVol][2];

    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                cross_array[i * GridLen * GridLen + j * GridLen + k][0] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]+ 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
                cross_array[i * GridLen * GridLen + j * GridLen + k][1] = 
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0]-
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * sc2[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1];
            }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k){
                cross_array[i * GridLen * GridLen + j * GridLen + k][0] = 
                cross_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k][0];
                cross_array[i * GridLen * GridLen + j * GridLen + k][1] = 
                cross_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k][1];
            }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k){
                cross_array[i * GridLen * GridLen + j * GridLen + k][0] *= PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k];
                cross_array[i * GridLen * GridLen + j * GridLen + k][1] *= PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k];
            }
    int klen = GridLen*sqrt(3.);
    int nk[klen];
    auto ccf = new double[klen];
    for(int i = 0; i < klen; ++i)
    {   
        nk[i] = 0;
        ccf[i] = 0;
    }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k){
                int kM = sqrt(i * i + j * j + k * k);
                ccf[kM] += sqrt(pow(cross_array[i * GridLen * GridLen + j * GridLen + k][0],2)+ 
                                pow(cross_array[i * GridLen * GridLen + j * GridLen + k][1],2));
                nk[kM] += 1;
            }
    delete[] cross_array;
    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        ccf[i] /= nk[i];
        ccf[i] /= pow(GridVol,2);
    }
    
    return ccf;
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
    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction (+:sum)
    #endif
    
    for(size_t i = 0; i < N; ++i)
    {
        sum += v0[i] * v1[i];
    }
    return sum;
}

// retrun sum of all elements of array w, w has length N
double array_sum(double* w, size_t N)
{
    double sum {0};
    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction (+:sum)
    #endif
    for(size_t i = 0; i < N; ++i)
    {
        sum += w[i];
    }
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

// given a Cloud-in-Cell interpreted 3d density array in fourier space, and the Gaussian kernel also in fourier space
// calculate the soothed density fileds and the Hessian elements of the gravitation potential fileds (df, xx, xy, xz, yy, yz, zz)
// return as pointer of <double> pointer cxx, with cxx[0] the address of df array, cxx[1] address of xx array and so on.
double** tidal_tensor(fftw_complex* sc, double* w)
{
    auto sc0 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    double** cxx = new double*[7];
    for(int i = 0; i < 7; ++i){
        cxx[i] = new double[GridVol];
    }
    fftw_plan pldf = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc0, cxx[0], FFTW_MEASURE);
    fftw_plan plxx = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[1], FFTW_MEASURE);
    fftw_plan plxy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[2], FFTW_MEASURE);
    fftw_plan plxz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[3], FFTW_MEASURE);
    fftw_plan plyy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[4], FFTW_MEASURE);
    fftw_plan plyz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[5], FFTW_MEASURE);
    fftw_plan plzz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[6], FFTW_MEASURE);

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    {
        sc0[i][0] = w[i] * sc[i][0];
        sc0[i][1] = w[i] * sc[i][1]; 
    }
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * i / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * i / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxx);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * j / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * j * j / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plyy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * j * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plyz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * k * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * k * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plzz);
    fftw_execute(pldf);
    
    fftw_destroy_plan(pldf);
    fftw_destroy_plan(plxx);
    fftw_destroy_plan(plxy);
    fftw_destroy_plan(plxz);
    fftw_destroy_plan(plyy);
    fftw_destroy_plan(plyz);
    fftw_destroy_plan(plzz);
    fftw_free(sc0);
    fftw_free(sc1);

    return cxx;
}

// solving a 3x3 real symmetric matrix and return the number of eigenvalue above a threshold Lambda_th=0
int eigen_classify(double xx, double xy, double xz, double yy, double yz, double zz)
{
    int n = 0;
    double lambda_th = 0;
    double t = xx + yy + zz;
    double mid1 = 2 * xx - yy - zz;
    double mid2 = 2 * yy - xx - zz;
    double mid3 = 2 * zz - xx - yy;
    double a = xx * xx + yy * yy + zz * zz - xx * yy - xx * zz - yy * zz + 3 * (xy * xy + xz * xz + yz * yz);
    double b = 9 * (mid1 * yz * yz + mid2 * xz * xz + mid3 * xy * xy) - 54 * xy * xz * yz - mid1 * mid2 * mid3;
    double phi = M_PI / 6;
    if(b > 0)
        phi = atan(sqrt(4 * a * a * a - b * b) / b) / 3.;
    else if(b < 0)
        phi = (atan(sqrt(4 * a * a * a - b * b) / b) + M_PI) / 3.;
    double lambda[3];
    lambda[0] = (t - 2 * sqrt(a) * cos(phi)) / 3;
    lambda[1] = (t - 2 * sqrt(a) * cos(phi + 2 * M_PI / 3)) / 3;
    lambda[2] = (t - 2 * sqrt(a) * cos(phi - 2 * M_PI / 3)) / 3;
    for(int i = 0; i < 3; ++i)
        if(lambda[i] > lambda_th)
            ++n;
    return n;
}

// from Hessian to web structure: 0-Voids, 1-sheets, 2-filaments, 3-knots
std::vector<int> web_classify(double** cxx, std::vector<Particle>& p)
{
    std::vector<int> s(p.size());
    #pragma omp parallel for
    for(int i = 0; i < p.size(); ++i)
    {
        int xs,ys,zs;   // BSpline have support [0,n+1],not centre in origin
        int64_t x,y,z,l;   
        xs = p[i].x / SimBoxL * GridLen + 0.5;
        ys = p[i].y / SimBoxL * GridLen + 0.5;
        zs = p[i].z / SimBoxL * GridLen + 0.5;
        x = (xs - 1) & (GridLen - 1);   // shift -1 for CIC, which is corresponding to BSpline n=1
        y = (ys - 1) & (GridLen - 1);
        z = (zs - 1) & (GridLen - 1);
        l = x * GridLen * GridLen + y * GridLen + z;
        s[i] = eigen_classify(cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], cxx[6][l]);
    }
    return s;
}

// from Hessian to web structure on grid: 0-Voids, 1-sheets, 2-filaments, 3-knots
std::vector<int> web_classify_to_grid(double** cxx)
{
    std::vector<int> s(GridVol);
    #pragma omp parallel for
    for(int64_t i = 0; i < GridLen; ++i)
        for(int64_t j = 0; j < GridLen; ++j)
            for(int64_t k = 0; k < GridLen; ++k)
            {
                int64_t x,y,z,l;   // BSpline have support [0,n+1],not centre in origin
                x = (i - 1) & (GridLen - 1);   // shift -1 for CIC, which is corresponding to BSpline n=1
                y = (j - 1) & (GridLen - 1);
                z = (k - 1) & (GridLen - 1);
                l = x * GridLen * GridLen + y * GridLen + z;
                s[i * GridLen * GridLen + j * GridLen + k] = eigen_classify(cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], cxx[6][l]);
            }
    return s;
}


// given dark matter fields dm and Gaussian smoothing radius Rs, return the web classify result at point p0
std::vector<int> environment(std::vector<Particle>& dm, double Rs, std::vector<Particle>& p0)
{
    force_base_type(0,1);
    force_kernel_type(2);
    auto begin = std::chrono::steady_clock::now();

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(2, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p0);
    
    fftw_free(sc);
    delete[] w;

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time difference EnvironmentClassify  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
    << "[ms]" << std::endl;

    return env;
}



double gaussian_radius_from_mass(double m_smooth) 
{
    return 1. / sqrt(TWOPI) * pow(m_smooth, 1./3); // rho_bar is needed: [m_smooth / rho_ar]
}

// p[i] = s1[i] * Hermitian[s2[i]]
fftw_complex* hermitian_product(fftw_complex* sc1, fftw_complex* sc2)
{
    auto c = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    {
        c[i][0] = sc1[i][0] * sc2[i][0] + sc1[i][1] * sc2[i][1];
        c[i][1] = sc1[i][1] * sc2[i][0] - sc1[i][0] * sc2[i][1];
    }

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

//=======================================================================================
//|||||||||||||||||||||||||||||||||||| Result interpret |||||||||||||||||||||||||||||||||
//=======================================================================================
void result_interpret(const double* s, std::vector<Particle>& p0, std::vector<double>& result)
{
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates

    auto begin4 = std::chrono::steady_clock::now();

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;


    for(size_t n = 0; n < p0.size(); ++n)
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
                    sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                            * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }
        result.push_back(sum);
    }

    auto end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 interpret = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;
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
double* count_in_sphere(const double R, std::vector<Particle>& p, std::vector<Particle>& p0)
{
    std::vector<Index> inner_index;
    std::vector<Index> cross_index;
    std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();

    fill_index_set(R/SimBoxL, inner_index, cross_index);

    auto count = new double[p0.size()];
    double temp;

    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction(+:temp)
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

    std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 count    = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return count;
}


// safe band is needed, i.e, no periodic boundary condictions added in
double* count_in_cylinder(double R, double H, std::vector<Particle>& p, std::vector<Particle>& p0)
{
    std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();
    auto count = new double[p0.size()];
    double temp;

    #pragma omp parallel for reduction(+:temp)
    for(size_t n = 0; n < p0.size(); ++n)
    {
        temp = 0;
        for(size_t i = 0; i < p.size(); ++i)
        {
            double xx = p[i].x - p0[n].x;
            double yy = p[i].y - p0[n].y;
            double zz = p[i].z - p0[n].z;
            if((abs(xx) < R) && (abs(yy) < R) && (abs(zz) < H/2))
                if(xx*xx + yy*yy  < R*R)
                    ++temp;
        }
            
        count[n]= temp;
    }

    std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 count_cyl = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return count;
}




// enviriamental parameter array locate in gride point, defined as the number of positive eigenvalue
// of tidle tensor of matter density fileds, which is obtand by Cloud-in-Cell interpolation of particles
// to grid point and then smoothed by a Gaussian kernel with radius R. The main process is working on 
// fourier space so we can take advantage of FFT, for detials see Hahn O., Porciani C., Carollo C. M., Dekel A., 2007, MNRAS, 375, 489
// https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..489H

// generate n random particle locate in box (0,boxsize)^3 
std::vector<Particle> default_random_particle(double boxsize, size_t n)
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> v(0,SimBoxL);

    std::vector<Particle> p0;
    for(size_t i = 0; i < n; ++i) 
        p0.push_back({v(e), v(e), v(e), 1.});
    
    return p0;
}

// x is the number of particles per side, i.e there are x^3 particles in all
// L is the box size and w is the safe band width, i.e random point locate in (w,L-w)^3
std::vector<Particle> generate_random_particle(int x, double L, double w)
{
    double diff;
    clock_t begin, end;

    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int64_t nps = x;
    const int64_t NPt = nps * nps * nps;
    const double boxL = L;
    const double safeband = w;

    begin = clock();
    for(int i = 0; i < nps; ++i)
        for(int j = 0; j < nps; ++j)
            for(int k = 0; k < nps; ++k)
                p.push_back({safeband + (i + u(e)) * (boxL - 2*safeband) / nps,
                              safeband + (j + u(e)) * (boxL - 2*safeband) / nps,
                              safeband + (k + u(e)) * (boxL - 2*safeband) / nps});
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "generating time for " << p.size() << " random points:  " << diff << "s" << std::endl;
    
    return p;
}

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

// nbin is the number of catalogue splited, return as vector of catalog
std::vector<std::vector<Particle>> halo_mass_split(std::vector<Particle>& hl, int nbin)
{
    std::vector<double> vecmass(hl.size());
    #pragma omp parallel for
    for(size_t i = 0; i < hl.size(); ++i) vecmass[i] = hl[i].weight;
    auto node = proto_sort(vecmass, nbin);

    std::vector<std::vector<Particle>> cata(nbin);
    for(auto x : hl){
        cata[classify_index(node,x.weight)].push_back(x);
    }
    for(auto x : cata) if(x.size() - hl.size()/nbin > nbin) {
        std::cout << "[func: SPLIT] !Warning, some nodes include multiple identical items\n";
        break;
    }
    for(auto x : cata) print_min_max_and_size(x);

    return cata;
}

void print_min_max_and_size(std::vector<Particle>& hl){
    double min{hl[0].weight},max{hl[0].weight};
    for(auto x : hl){
        if(x.weight > max) max = x.weight;
        else if(x.weight < min) min = x.weight;
    }
    std::cout << "size: " << hl.size() << ", min: " << min << ", max: " << max << std::endl; 
}
// which vector should trial been push back
int classify_index(std::vector<double>& node, double trial){
    int index{0};
    for(int i = 0; i < node.size(); ++i){
        if(trial > node[i]) ++index;
        else break;
    }
    return index;
}


// **************************************************************************
// return the nbin fraction node points of a double vector in ascending order
// **************************************************************************
std::vector<double> proto_sort(std::vector<double>& vec, int nbin)
{
    std::vector<double> node;
    auto max_id = maximum_index(vec);
    auto min_id = minimum_index(vec);

    double MAX{vec[max_id]}, MIN{vec[min_id]};
    double DELTA{MAX - MIN};
    const int REFINE{100};
    const size_t LINSIZE{vec.size() < 1e9 ? vec.size() * REFINE : vec.size()};
    const size_t VECLEN{vec.size()/nbin};

    auto count = new int[LINSIZE + 1](0);
    for(size_t i = 0; i < vec.size(); ++i){
        size_t idx = (vec[i] - MIN) / DELTA * LINSIZE;
        ++count[idx];
    }
    
    size_t node_id{0};
    for(int n = 0; n < nbin - 1; ++n){
        size_t sum{0}, finer{0};
        for(size_t i = node_id; i < LINSIZE + 1; ++i){
            sum += count[i];
            if(sum >= VECLEN) 
            {   
                node_id = i;
                break;
            }
        }
        if(count[node_id] > 1){
            finer =  VECLEN - (sum - count[node_id]);
            std::vector<double> nodex;
            for(size_t i = 0; i < vec.size(); ++i){
                size_t idx = (vec[i] - MIN) / DELTA * LINSIZE;
                if(idx == node_id) nodex.push_back(vec[i]);
            }
            auto index = limited_sort(nodex);
            node.push_back(nodex[index[finer-1]]);
            std::vector<double>().swap(nodex);
            std::vector<size_t>().swap(index);
            count[node_id] -= finer;
        }
        else 
            node.push_back(node_id * DELTA / LINSIZE + MIN);
        
    }
    delete count;

    return node;
}

// -------------------------------------------------------
// direct sort of double vector, return as ascending index
// -------------------------------------------------------
std::vector<size_t> limited_sort(std::vector<double> vec)
{
    std::cout << "[func: limited_sort] node size: " << vec.size() << std::endl;
    std::vector<size_t> sortedID;
    size_t max_id = maximum_index(vec);
    const double MAX{vec[max_id]};
    for(size_t i = 0; i < vec.size() - 1; ++i){
        size_t id = minimum_index(vec);
        vec[id] = MAX;
        sortedID.push_back(id);
    }
    sortedID.push_back(max_id);

    return sortedID;
}

size_t minimum_index(std::vector<double>& v)
{
    size_t idx{0};
    for(size_t i = 0; i < v.size(); ++i){
        if(v[i] < v[idx]) idx = i;
    }
    return idx;
}

size_t maximum_index(std::vector<double>& v)
{
    size_t idx{0};
    for(size_t i = 0; i < v.size(); ++i){
        if(v[i] > v[idx]) idx = i;
    }
    return idx;
}



std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl)
{
    auto dm_sc = sfc_r2c(sfc(dm),true);
    auto hl_sc = sfc_r2c(sfc(hl),true);

    double mh[GridLen](0);
    double mm[GridLen](0);
    double hh[GridLen](0);
    
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < (GridLen/2 + 1); ++k){
                auto ii = i < GridLen/2 + 1 ? i : GridLen - i;
                auto jj = j < GridLen/2 + 1 ? j : GridLen - j;
                int ll = sqrt(ii * ii + jj * jj + k * k);
                auto l = i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k;

                mh[ll] += dm_sc[l][0] * hl_sc[l][0] + dm_sc[l][1] * hl_sc[l][1];
                mm[ll] += dm_sc[l][0] * dm_sc[l][0] + dm_sc[l][1] * dm_sc[l][1];
                hh[ll] += hl_sc[l][0] * hl_sc[l][0] + hl_sc[l][1] * hl_sc[l][1];
            }
    std::vector<double> cross(GridLen/2 + 1);
    for(int i = 0; i < GridLen/2 + 1; ++i) cross[i] = mh[i]/sqrt(mm[i] * hh[i]);
    fftw_free(dm_sc);
    fftw_free(hl_sc);

    return cross;
}




/*
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

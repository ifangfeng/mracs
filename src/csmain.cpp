#include"csmain.h"



#define LOWER_RESOLUTION 8

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
    GridNum = 1UL << Resolution*3;

    RESOL = "L" + std::to_string(GridLen);
    RADII = "R" + std::to_string(Radius);
    GENUS = "DaubG" + std::to_string(phiGenus);
    
    if(BaseType == 0)
    {
        phi = B_Spline(phiGenus, SampRate);
        phi.push_back(0);
        phi.push_back(phiGenus + 1);
    }
    else if (BaseType == 1)
    {
        phi = Daubechies_Phi(phiGenus);
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

    //std::cout << "---Daubechies phi genus: " << phiGenus << std::endl;
    //std::cout << "---Wavelet phi supports: [0," << phiSupport << ") \n";
    //std::cout << "---Sampling points: " << phi.size() - 2 << std::endl;

    return phi;
}


//=======================================================================================
//||||||||||||||| sfc3d (scaling function coefficients of density field) ||||||||||||||||
//=======================================================================================
double* sfc_offset(std::vector<Particle>& p, Offset v)
{   
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates
    
    double total {0};
    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction (+:total)
    #endif
    for(auto x : p) total += x.weight;
    total /= p.size();

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
    {
        step[i] = i * SampRate;
    }

    auto s = new double[GridNum]();         // density field coefficients in v_j space
    

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    double s_temp[phiSupport * phiSupport * phiSupport];

    #ifdef IN_PARALLEL
    #pragma omp parallel for reduction (+:s_temp)
    #endif
    
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
                    phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]] * p[n].weight / total;
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
    size_t N_x= v.size()-2;
    double x0 = v[v.size()-2];
    double x1 = v[v.size()-1];
    

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
    size_t N_x = v.size()-2;
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
// particles first assigned to grid using different window function, then we can take
// advantage of FFT to get its Fourier Coefficiencs and average all orientation
// to have P(k) as function of scalar module k, notice that sfc function do the 
// assignment step exactly
//=======================================================================================
double* densityPowerFFT(double* s)
{
    auto sc = sfc_r2c(s);
    double* Pk_array = new double[GridNum];
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2);
            }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                Pk_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k];
            }
    fftw_free(sc);
    
    int klen = GridLen*sqrt(3.);
    int nk[klen];
    double* Pk = new double[klen];
    for(int i = 0; i < klen; ++i)
    {   
        nk[i] = 0;
        Pk[i] = 0;
    }
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k)
            {
                int kM = sqrt(i * i + j * j + k * k);
                Pk[kM] += Pk_array[i * GridLen * GridLen +j * GridLen + k];
                nk[kM] += 1;
            }
    delete[] Pk_array;
    
    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(GridNum,2);
    }
    
    return Pk;
}


//=======================================================================================
// calculate density power spectrum using projected density fileds mathematically
//=======================================================================================
double* densityPowerDWT(double* s)
{
    auto sc = sfc_r2c(s);
    double* Pk_array = new double[GridNum];
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0], 2) + 
                pow(sc[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1], 2);
            }
    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = GridLen/2 + 1; k < GridLen; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] = 
                Pk_array[((GridLen - i)%GridLen) * GridLen * GridLen + ((GridLen - j)%GridLen) * GridLen + GridLen - k];
            }
    fftw_free(sc);

    #ifdef IN_PARALLEL
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k)
            {
                Pk_array[i * GridLen * GridLen + j * GridLen + k] *=
                PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k];
            }
    int klen = GridLen*sqrt(3.);
    int nk[klen];
    double* Pk = new double[klen];
    for(int i = 0; i < klen; ++i)
    {   
        nk[i] = 0;
        Pk[i] = 0;
    }
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen; ++k)
            {
                int kM = sqrt(i * i + j * j + k * k);
                Pk[kM] += Pk_array[i * GridLen * GridLen +j * GridLen + k];
                nk[kM] += 1;
            }
    delete[] Pk_array;
    for(int i = 0; i < klen; ++i)
    {
        if(nk[i] != 0)
        Pk[i] /= nk[i];
        Pk[i] /= pow(GridNum,2);
    }
    
    return Pk;
}

//=======================================================================================
// calculate cross-correlation function c(k) of two density fileds in fourier space,
// using mass assignment and FFT method. where k is a scalar.
//=======================================================================================
double* densityCorrelationFFT(fftw_complex* sc1, fftw_complex* sc2)
{
    auto cross_array = new double[GridNum][2];

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
        ccf[i] /= pow(GridNum,2);
    }
    
    return ccf;
}


//=======================================================================================
// calculate cross-correlation function c(k) of two density fileds in fourier space,
// using projected density fileds mathematically, where k is a scalar.
//=======================================================================================
double* densityCorrelationDWT(fftw_complex* sc1, fftw_complex* sc2)
{
    auto cross_array = new double[GridNum][2];

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
        ccf[i] /= pow(GridNum,2);
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
        GridNum = 1UL << Resolution*3;
        delete[] PowerPhi;
        PowerPhi = PowerPhiFunc(GridLen);
        std::cout << "!MRACS resolution has been forced to " << j << "\n";
    }
}

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

// 'a' specify BaseType while 'n' phiGenus
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
            phi.push_back(0);
            phi.push_back(phiGenus + 1);
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

//=======================================================================================
//||||||||||||||| wfc3d (scaling function coefficients of window function) ||||||||||||||
//=======================================================================================
double* wfc(const double Radius, const double theta)
{
    const double DeltaXi = 1./GridLen;
    const double rescaleR {Radius * GridLen/SimBoxL};

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    auto WindowArray = new double[(GridLen+1) * (GridLen+1) * (GridLen+1)];
    auto w = new double[GridLen * GridLen * (GridLen/2+1)]();                             

    if(KernelFunc <= 2)
    {
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
        #ifdef IN_PARALLEL
        #pragma omp parallel for
        #endif
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                {
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = 
                    WindowFunction(rescaleR, i * DeltaXi, j * DeltaXi, k * DeltaXi) * PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k];
                }
    }
    else if(KernelFunc == 3)
    {
        double fz[GridLen+1];
        double dXitwo{pow(DeltaXi,2)};
        for(size_t i = 0; i <= GridLen; ++i) fz[i] = cos(TWOPI * rescaleR * cos(theta) * i*DeltaXi);
        double* fxy = new double[(GridLen+1) * (GridLen+1)];
        #ifdef IN_PARALLEL
        #pragma omp parallel for 
        #endif
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
            {
                fxy[i * (GridLen+1) + j] = std::cyl_bessel_j(0,TWOPI*sin(theta)*rescaleR*sqrt(i*i*dXitwo+j*j*dXitwo));
            }
        #ifdef IN_PARALLEL
        #pragma omp parallel for 
        #endif
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                {
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] = fxy[i * (GridLen+1) + j] * fz[k]; 
                    //WindowFunction_Dual_Ring(rescaleR, theta, i * DeltaXi, j * DeltaXi, k * DeltaXi);
                }
        #ifdef IN_PARALLEL
        #pragma omp parallel for 
        #endif
        for(size_t i = 0; i <= GridLen; ++i)
            for(size_t j = 0; j <= GridLen; ++j)
                for(size_t k = 0; k <= GridLen; ++k)
                {
                    WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k] *= PowerPhi[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k];
                }
        delete[] fxy;
    }
    #ifdef IN_PARALLEL     
    #pragma omp parallel for
    #endif
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2+1; ++k)
            {
                w[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k]
                = WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k]
                + WindowArray[(GridLen-i) * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + k]
                + WindowArray[i * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + k]
                + WindowArray[i * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + (GridLen-k)]
                + WindowArray[(GridLen-i) * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + k]
                + WindowArray[(GridLen-i) * (GridLen+1) * (GridLen+1) + j * (GridLen+1) + (GridLen-k)]
                + WindowArray[i * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + (GridLen-k)]
                + WindowArray[(GridLen-i) * (GridLen+1) * (GridLen+1) + (GridLen-j) * (GridLen+1) + (GridLen-k)];
            }

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d    = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()
    << "[ms]" << std::endl;

    delete[] WindowArray;

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

fftw_complex* sfc_r2c(double* s)
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

    return sc;
}

double* convol_c2r(fftw_complex* sc, double* w)
{
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto c = new double[GridNum];

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
    for(size_t i = 0; i < GridNum; ++i) c[i] /= GridNum;

    fftw_free(sc1);
    return c;
}

char* tidal_tensor(fftw_complex* sc, double* w)
{
    auto sc0 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    double* cxx[3][3] = {nullptr};
    for(int xi = 0; xi < 3; ++xi)
        for(int xj =0; xj < 3; ++xj){
            cxx[xi][xj] = new double[GridNum];
        }
    fftw_plan plxx = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[0][0], FFTW_MEASURE);
    fftw_plan plxy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[0][1], FFTW_MEASURE);
    fftw_plan plxz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[0][2], FFTW_MEASURE);
    fftw_plan plyy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[1][1], FFTW_MEASURE);
    fftw_plan plyz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[1][2], FFTW_MEASURE);
    fftw_plan plzz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[2][2], FFTW_MEASURE);

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
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = i * i / (i*i + j*j + k*k); 
            }
    fftw_execute(plxx);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = i * j / (i*i + j*j + k*k); 
            }
    fftw_execute(plxy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = i * k / (i*i + j*j + k*k); 
            }
    fftw_execute(plxz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = j * j / (i*i + j*j + k*k); 
            }
    fftw_execute(plyy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = j * k / (i*i + j*j + k*k); 
            }
    fftw_execute(plyz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * k * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = k * k / (i*i + j*j + k*k); 
            }
    fftw_execute(plzz);
    
    auto c = new char[GridNum];

    for(size_t i = 0; i < GridNum; ++i)
    {

    }
    
    return nullptr;
}


// void solving_and_classifying(double* a[3][3], char* cptr, size_t N)
// {
//     for(size_t i = 0; i < N; ++i)
//     {
//         double lambda0,lambda1,lambda2;
//         a[0][0][i],a[0][0][i],a[0][0][i],a[0][0][i],a[0][0][i],a[0][0][i],
//     }
// }
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
        phi = atan(sqrt(4*a*a*a - b*b) / b) / 3.;
    else if(b < 0)
        phi = (atan(sqrt(4*a*a*a - b*b) / b) + M_PI) / 3.;
    double lambda[3];
    lambda[0] = (t - 2*sqrt(a)*cos(phi)) / 3;
    lambda[1] = (t - 2*sqrt(a)*cos(phi + 2*M_PI / 3)) / 3;
    lambda[1] = (t - 2*sqrt(a)*cos(phi - 2*M_PI / 3)) / 3;
    for(int i = 0; i < 3; ++i)
        if(lambda[i] > lambda_th)
            ++n;
    return n;
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
double* convol3d(double* s, double* w)
{
    auto begin3 = std::chrono::steady_clock::now();

    auto sc = sfc_r2c(s);
    auto c = convol_c2r(sc, w);
    auto end3 = std::chrono::steady_clock::now();

    std::cout << "Time difference 3 convl3d  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count()
    << "[ms]" << std::endl;

    fftw_free(sc);
    return c;
}

//=======================================================================================
//|||||||||||||||||||||||||||||||||||| Result interpret |||||||||||||||||||||||||||||||||
//=======================================================================================
void result_interpret(const double* s, std::vector<Particle>& p0, std::vector<double>& result)
{
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const int SampRate   = (phi.size() - 2) / phiSupport;
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


double* project_value(const double* s, std::vector<Particle>& p0)
{
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const double ScaleFactor {GridLen/SimBoxL};   //used to rescale particle coordinates

    auto begin4 = std::chrono::steady_clock::now();

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

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
    std::cout << "Time difference 4 interpret = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

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

// enviriamental parameter array locate in gride point, defined as the number of positive eigenvalue
// of tidle tensor of matter density fileds, which is obtand by Cloud-in-Cell interpolation of particles
// to grid point and then smoothed by a Gaussian kernel with radius R. The main process is working on 
// fourier space so we can take advantage of FFT, for detials see Hahn O., Porciani C., Carollo C. M., Dekel A., 2007, MNRAS, 375, 489
// https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..489H










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

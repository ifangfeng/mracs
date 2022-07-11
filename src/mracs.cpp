#include"mracs.hpp"


int main()
{

    welcome();

    read_parameter();

    auto phi = read_in_phi(phiGenus, DIREC);
    auto p   = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    auto s   = scaling_function_coefficients(p, phi, Resolution, SimBoxL);
    auto w   = window_function_coefficients(phi, Resolution, SimBoxL, Radius);

    specialized_convolution_3d(s, w, Resolution);
    
    result_interpret(s, phi, Resolution, SimBoxL, p.size());

    return 0;
}





    /*
    while(iprmfs >> itemp)
    {
        if(itemp=="Resolution")
        {
            iprmfs >> temp;
            Resolution = atoi(temp.c_str());
        }
        else if(itemp=="phiGenus")
        {
            iprmfs >> temp;
            phiGenus = atoi(temp.c_str());
        }
        else if(itemp=="SampRate")
        {
            iprmfs >> temp;
            SampRate = atoi(temp.c_str());
        }
        else if(itemp == "KernelFunc")
        {
            iprmfs >> temp;
            KernelFunc = atoi(temp.c_str());
        }
        else if(itemp == "Radius")
        {
            iprmfs >> temp;
            double Radius = atof(temp.c_str());
        }
        else if(itemp == "DIREC")
        {
            iprmfs >> temp;
            std::string DIREC = temp;
        }
        else if(itemp == "MilleCata")
        {
            iprmfs >> temp;
            std::string MilleCata = temp;
        }
        else if(itemp == "SimBoxL")
        {
            iprmfs >> temp;
            double SimBoxL  = atof(temp.c_str());
        }
        
    

    std::cout << "Reading para.txt " << std::endl;
    std::cout << "-> Resolution =  " << Resolution << std::endl;
    std::cout << "-> phiGenus   =  " << phiGenus << std::endl;
    std::cout << "-> SampRate   =  " << SampRate << std::endl;
    std::cout << "-> KernelFunc =  " << KernelFunc << std::endl;
    std::cout << "-> Radius     =  " << Radius << std::endl;
    std::cout << "-> DIREC      =  " << DIREC << std::endl;
    std::cout << "-> MilleCata  =  " << MilleCata << std::endl;
    std::cout << "-> SimBoxL    =  " << SimBoxL << std::endl;

    */



/*

// in case the .bin file store in diffrent endianness
// it happens Millennium Run galaxy catalog use Big Endian
template<class T> void readBigEndian(std::ifstream& i, T& a)
{
    for(int n = sizeof(a) - 1; n >= 0 ; --n)
    {
        i.read(((char*) &a) + n, sizeof(char));
    }
}



//=======================================================================================
//|||||||||||||||||||||||||||||||| read in particle data ||||||||||||||||||||||||||||||||
//=======================================================================================
//----- Millennium Run galaxy catalog, croton_etal
//----- https://wwwmpa.mpa-garching.mpg.de/galform/agnpaper/
//---------------------------------------------------------------------------------------
std::vector<Galaxy> read_in_Millennium_Run_galaxy_catalog(const std::string MilleCata)
{
    std::string ipfname {MilleCata};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs) 
    {   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }


    std::vector<Galaxy> p;

    std::chrono::steady_clock::time_point begin0, end0;

    std::string extension = ipfname.substr(ipfname.find_last_of(".") + 1);
    std::string ministr   = ipfname.substr(ipfname.find_last_of(".") - std::strlen("mini"), std::strlen("mini"));

    int miniflag = (ministr=="mini") ? 1 : 0;
    
    //---------------------------------------------------------------------------------------
    //================================== read if in ASCII ===================================
    //---------------------------------------------------------------------------------------
    if(extension=="ascii"||extension=="ASCII")
    {
        std::cout << "Reading ASCII file..." << std::endl;
        begin0 = std::chrono::steady_clock::now();

        float x, y, z;                              // position in Mpc/h;
        float vx, vy, vz;                           // velocity in km/s;
        float Mag_u, Mag_g,                         // total galaxy magnitudes in ----
            Mag_r, Mag_i, Mag_z;                    // ---- (AB) standard SDSS filters
        float BulgeMag_u, BulgeMag_g,
            BulgeMag_r, BulgeMag_i, BulgeMag_z;     // bulge magnitude only
        float StellarMass, BulgeMass, ColdGas, 
            HotGas, EjectedMass, BlackHoleMass, Sfr;// all mass in 10^10Msun/h, Sfr in Msun/yr


        while(ifs >> x >> y >> z >> vx >> vy >> vz >> Mag_u >> Mag_g >> Mag_r >> Mag_i >> Mag_z >>
            BulgeMag_u >> BulgeMag_g >> BulgeMag_r >> BulgeMag_i >> BulgeMag_z >> StellarMass >>
            BulgeMass >> ColdGas >> HotGas >> EjectedMass >> BlackHoleMass >> Sfr)
        {
            p.push_back(Galaxy{x, y, z, vx, vy, vz, Mag_u, Mag_g, Mag_r, Mag_i, Mag_z, BulgeMag_u,
                                BulgeMag_g, BulgeMag_r, BulgeMag_i, BulgeMag_z, StellarMass,
                                BulgeMass, ColdGas, HotGas, EjectedMass, BlackHoleMass, Sfr});
        }

        end0 = std::chrono::steady_clock::now();
    }
    //---------------------------------------------------------------------------------------
    //===================================== if in binary ====================================
    //---------------------------------------------------------------------------------------
    else if(extension=="bin"||extension=="BIN")
    {
        std::cout << "Reading binary file..." << std::endl;
        begin0 = std::chrono::steady_clock::now();
        
        unsigned int Total;                         // read in header of file
        if(miniflag)
        {
            unsigned short x;
            readBigEndian(ifs, x);
            Total = x;
        }
        else
        {
            readBigEndian(ifs, Total);
        }

        p.resize(Total);                            // vector of Particles p
        
        Galaxy x;
        int L = x.size();                           // number of parameter of each Galaxy
        float* a = new float[Total * L];

        for(int i = 0; i < L; ++i)
            for(int j = 0; j < Total; ++j)
                {
                    readBigEndian(ifs, a[i+j*L]);   // read all particle data
                }                                   // it takes seconds of time
        
        for(int i = 0; i < Total; ++i)
        {
            p[i].init(&a[i*L]);                     // assign to vector p
        }

        delete[] a;
        end0 = std::chrono::steady_clock::now();
    }
    //---------------------------------------------------------------------------------------
    //===================================== invalid file ====================================
    //---------------------------------------------------------------------------------------
    else
    {
        std::cout << "<" + ipfname + "> File has no valid extension, Abort " << std::endl;
        std::terminate();
    }


    std::cout << "---Number of Particles: " << p.size() << std::endl;

    std::cout << "Time difference 0 Read Particles = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}


// needed for binary I/O
template<class T> char* as_bytes(T& i)
{
    void* addr = &i;                    // get the address of the first byte
                                        // of memory used to store the object
    return static_cast<char*>(addr);    // treat that memory as bytes
}

//=======================================================================================
//|||||||||||||||||||||||||||||| read in wavelets phi data ||||||||||||||||||||||||||||||
//=======================================================================================
//---- numerical value of Daubechies scaling function of genus X, locate in closed
//---- interval [0, 2X-1] with sampling rate 1000 points per unit length (points / 1)
//---------------------------------------------------------------------------------------
std::vector<double> read_in_phi(const int phiGenus, const std::string DIREC)
{
    const int phiStart   {0};                     //Wavelet Phi0() has compact support, 
    const int phiEnd     {2*phiGenus - 1};        //start in x == 0, end in phi_end == 2n-1
    const int phiSupport {phiEnd - phiStart};     //Wavelet Phi0() has compact support 
    const int SampRate   {1000};                  //Wavelet Phi sampling rate (points / 1)

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;
    
    std::string iwfname {DIREC + "/wfile/DaubechiesG" + std::to_string(phiGenus) + "Phi.bin"};
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
template<class T> double* scaling_function_coefficients(std::vector<T>& p, std::vector<double>& phi, const int J , const double SimBoxL)
{
    const int L          = 1 << J;
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const int SampRate   = (phi.size() - 2) / phiSupport;
    const double ScaleFactor {L/SimBoxL};   //used to rescale particle coordinates


    auto s = new double[L*L*L]();         // density field coefficients in v_j space

    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    for(int n = 0; n < p.size(); ++n)
    {
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
                    s[((xxc - i) & (L - 1)) * L * L+ ((yyc - j) & (L - 1)) * L + ((zzc - k) & (L - 1))] 
                    += phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }
                    
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 1 sfc3d = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    return s;
}

//read in Real Value Vector v, return an array s, which store the continuous fourier
//transform of v, located in frequency space [k0, k1) with N_k sampling points
double* Spectrum1(std::vector<double>& v, double k0, double k1, int N_k)
{
    int N_x= v.size()-2;
    double x0 = v[v.size()-2];
    double x1 = v[v.size()-1];
    
    //const double TWOPI {2*M_PI};
    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};

    double Real, Image, Temp;
    double Phase = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    
    double* s = new double[N_k * 2]();
    for(int i = 0; i < N_k; ++i)
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
    
    //const double TWOPI {2*M_PI};
    const double Delta_x {(x1-x0)/N_x};
    const double Delta_k {(k1-k0)/N_k};

    double Real, Image, Temp, Phase;
    double PhaseOut = -TWOPI * k0 * x0;
    double DeltaPhase = -TWOPI * k0 * Delta_x;
    const double dPhaseOut = -TWOPI * Delta_k * x0;
    const double dDeltaPhase = -TWOPI * Delta_k * Delta_x;
    
    double* s = new double[N_k * 2]();
    for(int i = 0; i < N_k; ++i)
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
//frequency space [k0, k1) with N_k sampling points, return as array p[]
double* PowerSpectrum(std::vector<double>& v, double k0, double k1, int N_k)
{
    double* s = Spectrum(v, k0, k1, N_k);
    double* p = new double[N_k];
    for(int i = 0; i < N_k; ++i)
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
    const double rescaleR {Radius * L / SimBoxL};

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    double* PowerPhi = PowerSpectrum(phi, 0, 1, L);

    auto w = new double[L*L*L]();         //window function coefficients in v_j space

    const double DeltaXi = 1./(L);
    for(int i = 0; i < L; ++i)
        for(int j = 0; j < L; ++j)
            for(int k = 0; k < L; ++k)
                w[i * L * L + j * L + k] = PowerPhi[i] * PowerPhi[j] * PowerPhi[k]
                      * WindowFunction_Sphere(rescaleR, i*DeltaXi, j*DeltaXi, k*DeltaXi);
    

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()
    << "[ms]" << std::endl;

    return w;
}


//calculate two array's inner product and store result in the first one
template<class T> void inner_product0(T& v0, T& v1, int N)
{
    for(int i = 0; i < N; ++i)
        v0[i] *= v1[i];
}

//s and w are 3d real array, s in physical space while w in frequency space
//convol == fftback(inner_product(fft(s), w)), Matrix3D == L*L*L , L == 2^J
//convolution result store in s
void specialized_convolution_3d(double* s, double* w, int J)
{
    const int L {1<<J};
    const int N {L*L*L};

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();

    auto sc = new double[N*2];
    for(int i = 0; i < N; ++i)
        sc[2*i] = s[i];
    FFT3D_CUBIC(sc, J, 1);
    for(int i = 0; i < N; ++i)
    {
        sc[2*i] *= w[i];
        sc[2*i+1] *= w[i]; 
    }
    FFT3D_CUBIC(sc, J, 0);
    for(int i = 0; i < N; ++i)
        s[i] = sc[2*i]/N;

    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3 convol3d = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count()
    << "[ms]" << std::endl;

    delete[] sc;
}
*/












/*
//=======================================================================================
//|||||||||||||||||||||||||||||| read in wavelets phi data ||||||||||||||||||||||||||||||
//=======================================================================================
//---- numerical value of Daubechies scaling function of genus X, locate in closed
//---- interval [0, 2X-1] with sampling rate 1000 points per unit length (points / 1)
//---------------------------------------------------------------------------------------
vector<double> read_in_phi(const int phiGenus, const string DIREC)
{
    const int phiStart   {0};                     //Wavelet Phi0() has compact support, 
    const int phiEnd     {2*phiGenus - 1};        //start in x == 0, end in phi_end == 2n-1
    const int phiSupport {phiEnd - phiStart};     //Wavelet Phi0() has compact support 
    const int SampRate   {1000};                  //Wavelet Phi sampling rate (points / 1)

    vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;
    
    string iwfname {DIREC + "/wfile/DaubechiesG" + to_string(phiGenus) + "Phi.bin"};
    ifstream iwfs{iwfname, ios_base::binary};
    if(!iwfs) 
    {   
        cout << "!Reading " + iwfname + " with error..." << endl;
        terminate();
    }
    
    double temp;
    vector<double> phi;
    while (iwfs.read(as_bytes(temp), sizeof(double)))
    {
        phi.push_back(temp);
    }
    
    phi.push_back(phiStart);
    phi.push_back(phiEnd);

    std::cout << "---Daubechies phi genus: " << phiGenus << std::endl;
    std::cout << "---Wavelet phi supports: [0," << phiSupport << ") \n";
    std::cout << "---Sampling points: " << phi.size() - 2 << '\n';

    return phi;
}

//=======================================================================================
//|||||||||||||||||||||||||||||||| read in particle data ||||||||||||||||||||||||||||||||
//=======================================================================================
//----- Millennium Run galaxy catalog, croton_etal
//----- https://wwwmpa.mpa-garching.mpg.de/galform/agnpaper/
//---------------------------------------------------------------------------------------
vector<Galaxy> read_in_Millennium_Run_galaxy_catalog(const string MilleCata)
{
    string ipfname {MilleCata};
    ifstream ifs {ipfname, ios_base::binary};
    if(!ifs) 
    {   
        cout << "!Reading " + ipfname + " with error..." << endl;
        terminate();
    }


    vector<Galaxy> p;

    std::chrono::steady_clock::time_point begin0, end0;

    string extension = ipfname.substr(ipfname.find_last_of(".") + 1);
    string ministr   = ipfname.substr(ipfname.find_last_of(".") - strlen("mini"), strlen("mini"));

    int miniflag = (ministr=="mini") ? 1 : 0;
    
    //---------------------------------------------------------------------------------------
    //================================== read if in ASCII ===================================
    //---------------------------------------------------------------------------------------
    if(extension=="ascii"||extension=="ASCII")
    {
        cout << "Reading ASCII file..." << endl;
        begin0 = std::chrono::steady_clock::now();

        float x, y, z;                              // position in Mpc/h;
        float vx, vy, vz;                           // velocity in km/s;
        float Mag_u, Mag_g,                         // total galaxy magnitudes in ----
            Mag_r, Mag_i, Mag_z;                    // ---- (AB) standard SDSS filters
        float BulgeMag_u, BulgeMag_g,
            BulgeMag_r, BulgeMag_i, BulgeMag_z;     // bulge magnitude only
        float StellarMass, BulgeMass, ColdGas, 
            HotGas, EjectedMass, BlackHoleMass, Sfr;// all mass in 10^10Msun/h, Sfr in Msun/yr


        while(ifs >> x >> y >> z >> vx >> vy >> vz >> Mag_u >> Mag_g >> Mag_r >> Mag_i >> Mag_z >>
            BulgeMag_u >> BulgeMag_g >> BulgeMag_r >> BulgeMag_i >> BulgeMag_z >> StellarMass >>
            BulgeMass >> ColdGas >> HotGas >> EjectedMass >> BlackHoleMass >> Sfr)
        {
            p.push_back(Galaxy{x, y, z, vx, vy, vz, Mag_u, Mag_g, Mag_r, Mag_i, Mag_z, BulgeMag_u,
                                BulgeMag_g, BulgeMag_r, BulgeMag_i, BulgeMag_z, StellarMass,
                                BulgeMass, ColdGas, HotGas, EjectedMass, BlackHoleMass, Sfr});
        }

        end0 = std::chrono::steady_clock::now();
    }
    //---------------------------------------------------------------------------------------
    //===================================== if in binary ====================================
    //---------------------------------------------------------------------------------------
    else if(extension=="bin"||extension=="BIN")
    {
        cout << "Reading binary file..." << endl;
        begin0 = std::chrono::steady_clock::now();
        
        unsigned int Total;                         // read in header of file
        if(miniflag)
        {
            unsigned short x;
            readBigEndian(ifs, x);
            Total = x;
        }
        else
        {
            readBigEndian(ifs, Total);
        }

        p.resize(Total);                            // vector of Particles p
        
        Galaxy x;
        int L = x.size();                           // number of parameter of each Galaxy
        float* a = new float[Total * L];

        for(int i = 0; i < L; ++i)
            for(int j = 0; j < Total; ++j)
                {
                    readBigEndian(ifs, a[i+j*L]);   // read all particle data
                }                                   // it takes seconds of time
        
        for(int i = 0; i < Total; ++i)
        {
            p[i].init(&a[i*L]);                     // assign to vector p
        }

        delete[] a;
        end0 = std::chrono::steady_clock::now();
    }
    //---------------------------------------------------------------------------------------
    //===================================== invalid file ====================================
    //---------------------------------------------------------------------------------------
    else
    {
        cout << "<" + ipfname + "> File has no valid extension, Abort " << endl;
        terminate();
    }


    std::cout << "---Number of Particles: " << p.size() << std::endl;

    std::cout << "Time difference 0 Read Particles = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << endl;

    return p;
}


//=======================================================================================
//||||||||||||||| sfc3d (scaling function coefficients of density field) ||||||||||||||||
//=======================================================================================
template<class T> double* scaling_function_coefficients(vector<T>& p, vector<double>& phi, const int J , const double SimBoxL)
{
    const int L          = 1 << J;
    const int phiStart   = phi[phi.size() - 2];
    const int phiEnd     = phi[phi.size() - 1];
    const int phiSupport = phiEnd - phiStart;
    const int SampRate   = (phi.size() - 2) / phiSupport;
    const double ScaleFactor {L/SimBoxL};   //used to rescale particle coordinates


    auto s = new double[L*L*L]();         // density field coefficients in v_j space

    vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i)
        step[i] = i * SampRate;

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    for(int n = 0; n < p.size(); ++n)
    {
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
                    s[((xxc - i) & (L - 1)) * L * L+ ((yyc - j) & (L - 1)) * L + ((zzc - k) & (L - 1))] 
                    += phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
                }
                    
    }

    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 1 sfc3d = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << endl;

    return s;
}


//=======================================================================================
//||||||||||||||| wfc3d (scaling function coefficients of window function) ||||||||||||||
//=======================================================================================
double* window_function_coefficients(vector<double>& phi, const int J, const double SimBoxL, const double Radius)
{
    const int L {1<<J};
    const double rescaleR {Radius * L / SimBoxL};

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

    double* PowerPhi = PowerSpectrum(phi, 0, 1, L);

    auto w = new double[L*L*L]();         //window function coefficients in v_j space

    const double DeltaXi = 1./(L);
    for(int i = 0; i < L; ++i)
        for(int j = 0; j < L; ++j)
            for(int k = 0; k < L; ++k)
                w[i * L * L + j * L + k] = PowerPhi[i] * PowerPhi[j] * PowerPhi[k]
                      * WindowFunction_Sphere(rescaleR, i*DeltaXi, j*DeltaXi, k*DeltaXi);
    

    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time difference 2 wfc3d = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()
    << "[ms]" << endl;

    return w;
}


//calculate two array's inner product and store result in the first one
template<class T> void inner_product0(T& v0, T& v1, int N)
{
    for(int i = 0; i < N; ++i)
        v0[i] *= v1[i];
}

//s and w are 3d real array, s in physical space while w in frequency space
//convol == fftback(inner_product(fft(s), w)), Matrix3D == L*L*L , L == 2^J
//convolution result store in s
void specialized_convolution_3d(double* s, double* w, int J)
{
    const int L {1<<J};
    const int N {L*L*L};

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();

    auto sc = new double[N*2];
    for(int i = 0; i < N; ++i)
        sc[2*i] = s[i];
    FFT3D_CUBIC(sc, J, 1);
    for(int i = 0; i < N; ++i)
    {
        sc[2*i] *= w[i];
        sc[2*i+1] *= w[i]; 
    }
    FFT3D_CUBIC(sc, J, 0);
    for(int i = 0; i < N; ++i)
        s[i] = sc[2*i]/N;

    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3 convol3d = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count()
    << "[ms]" << endl;

    delete[] sc;
}
*/



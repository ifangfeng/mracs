#include"MRACS_Readin.h"



//=======================================================================================
//|||||||||||||||||||||||||||||||| read in particle data ||||||||||||||||||||||||||||||||
//=======================================================================================
//----- Millennium Run galaxy catalog, croton_etal
//----- https://wwwmpa.mpa-garching.mpg.de/galform/agnpaper/
//---------------------------------------------------------------------------------------
std::vector<Galaxy> read_in_Millennium_Run_galaxy_catalog(const std::string DataDirec)
{
    std::string ipfname {DataDirec};
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

    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}

// read in TNG dumped velocities or coordinates data, double precision
std::vector<Particle> read_in_TNG_3vector(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<Particle> p;
    double a[3];
    void* addr = a;

    while(ifs.read(static_cast<char*>(addr), 3*sizeof(double)))
        p.push_back({a[0],a[1],a[2],1.});


    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}



// float x, y, z
std::vector<Particle> read_in_DM_3vector(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<Particle> p;
    float a[3];
    void* addr = a;

    while(ifs.read(static_cast<char*>(addr), 3*sizeof(float)))
        p.push_back({a[0],a[1],a[2],1.});


    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}
// float x, y, z, M_virial 
std::vector<Particle> read_in_Halo_4vector(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<Particle> p;
    float a[4];
    void* addr = a;

    while(ifs.read(static_cast<char*>(addr), 4*sizeof(float)))
        p.push_back({a[0],a[1],a[2],a[3]});

    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}

// float x, y, z
std::vector<Particle> read_in_Halo_3vector(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<Particle> p;
    float a[4];
    void* addr = a;

    while(ifs.read(static_cast<char*>(addr), 4*sizeof(float)))
        p.push_back({a[0],a[1],a[2],1.});

    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}

// float x, y, z, Mass_vir, concentration, spin
std::vector<Halo> read_in_Halo_6vector(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<Halo> p;
    float a[6];
    void* addr = a;

    while(ifs.read(static_cast<char*>(addr), 6*sizeof(float)))
        p.push_back({a[0],a[1],a[2],a[3],a[4],a[5]});

    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}


//
std::vector<double> read_in_float(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<double> p;
    float a;
    void* addr = &a;

    while(ifs.read(static_cast<char*>(addr), sizeof(float)))
        p.push_back(a);

    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}

std::vector<double> read_in_double(std::string DataDirec){
    std::string ipfname {DataDirec};
    std::ifstream ifs {ipfname, std::ios_base::binary};
    if(!ifs){   
        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
        std::terminate();
    }
    std::chrono::steady_clock::time_point begin0, end0;
    std::cout << "Reading binary file..." << std::endl;
    begin0 = std::chrono::steady_clock::now();
    std::vector<double> p;
    double a;
    void* addr = &a;

    while(ifs.read(static_cast<char*>(addr), sizeof(double)))
        p.push_back(a);

    end0 = std::chrono::steady_clock::now();

    std::cout << "---Number of Particles: " << p.size() << std::endl;
    std::cout << "Time difference 0 Read_in  = " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
    << "[ms]" << std::endl;

    return p;
}

//template<class T> std::vector<T> read_in_1vector(std::string DataDirec){
//    std::string ipfname {DataDirec};
//    std::ifstream ifs {ipfname, std::ios_base::binary};
//    if(!ifs){   
//        std::cout << "!Reading " + ipfname + " with error..." << std::endl;
//        std::terminate();
//    }
//    std::chrono::steady_clock::time_point begin0, end0;
//    std::cout << "Reading binary file..." << std::endl;
//    begin0 = std::chrono::steady_clock::now();
//    std::vector<T> p;
//    T a;
//    void* addr = &a;
//
//    while(ifs.read(static_cast<char*>(addr), sizeof(T)))
//        p.push_back(a);
//
//    end0 = std::chrono::steady_clock::now();
//
//    std::cout << "---Number of Particles: " << p.size() << std::endl;
//    std::cout << "Time difference 0 Read_in  = " 
//    << std::chrono::duration_cast<std::chrono::milliseconds>(end0 - begin0).count()
//    << "[ms]" << std::endl;
//
//    return p;
//}
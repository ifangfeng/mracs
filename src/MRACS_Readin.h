#include"MRACS_Main.h"

#ifndef MRACS_READIN
#define MRACS_READIN

// needed for binary I/O
template<class T> char* as_bytes(T& i)
{
    void* addr = &i;                            // get the address of the first byte
                                                // of memory used to store the object
    return static_cast<char*>(addr);            // treat that memory as bytes
}

// in case the .bin file store in diffrent endianness
template<class T> void readBigEndian(std::ifstream& i, T& a)
{
    for(int n = sizeof(a) - 1; n >= 0 ; --n)
    {
        i.read(((char*) &a) + n, sizeof(char));
    }
}
// Millennium Run galaxy catalog
struct Galaxy
{
    float x, y, z;                              // position in Mpc/h;
    float vx, vy, vz;                           // velocity in km/s;
    float Mag_u, Mag_g, Mag_r, Mag_i, Mag_z;    // total galaxy magnitudes in (AB) standard SDSS filters
    float BulgeMag_u, BulgeMag_g,
           BulgeMag_r, BulgeMag_i, BulgeMag_z;  // bulge magnitude only
    float StellarMass, BulgeMass, ColdGas, 
           HotGas, EjectedMass, BlackHoleMass, Sfr; // all mass in 10^10Msun/h, Sfr in Msun/yr

    int size() {return 23;}                     // number of parameter {x,y,x,vx,vy...} == 23

    // initialize Galaxy from contineous memory <float*> a
    void init(float* a)                         
    {
        x               =   a[0];
        y               =   a[1];
        z               =   a[2];
        vx              =   a[3];
        vy              =   a[4];
        vz              =   a[5];
        Mag_u           =   a[6];
        Mag_g           =   a[7];
        Mag_r           =   a[8];
        Mag_i           =   a[9];
        Mag_z           =   a[10];
        BulgeMag_u      =   a[11];
        BulgeMag_g      =   a[12];
        BulgeMag_r      =   a[13];
        BulgeMag_i      =   a[14];
        BulgeMag_z      =   a[15];
        StellarMass     =   a[16];
        BulgeMass       =   a[17];
        ColdGas         =   a[18];
        HotGas          =   a[19];
        EjectedMass     =   a[20];
        BlackHoleMass   =   a[21];
        Sfr             =   a[22];
    }
};
std::vector<Galaxy> read_in_Millennium_Run_galaxy_catalog(const std::string DataDirec);
std::vector<Particle> read_in_TNG_3vector(std::string DataDirec);
std::vector<Particle> read_in_DM_3vector(std::string DataDirec);
std::vector<Particle> read_in_Halo_4vector(std::string DataDirec);
std::vector<Particle> read_in_Halo_3vector(std::string DataDirec);



#endif
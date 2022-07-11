#include"mracs.hpp"

//#define UNIFORM
#define GRIDSIZE 128
#define TESTPOINTS 100


int main()
{
    read_parameter();
    auto phi = read_in_phi(phiGenus, DIREC);
    
    double scalefactor = SimBoxL/GRIDSIZE;
    std::vector<Particle> p;

    #ifdef UNIFORM
    std::cout << "UNIFORM: " << std::endl;
    std::cout << "-Num of Particles: " << GRIDSIZE << "^3" << std::endl;
    for(int i = 0; i < GRIDSIZE; ++i)
        for(int j = 0; j < GRIDSIZE; ++j)
            for(int k = 0; k < GRIDSIZE; ++k)
              p.push_back(Particle{i*scalefactor, j*scalefactor, k*scalefactor, 1.});
    #else
    auto g   = read_in_Millennium_Run_galaxy_catalog(MilleCata);
    for(auto i : g) p.push_back(Particle{i.x, i.y, i.z, 1.});
    #endif

    auto s   = scaling_function_coefficients(p, phi, Resolution, SimBoxL);
    auto w   = window_function_coefficients(phi, Resolution, SimBoxL, Radius);
    specialized_convolution_3d(s, w, Resolution);

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);

    std::vector<Particle> p0;

    for(int i = 0; i < TESTPOINTS; ++i)
    {   
        p0.push_back(Particle{u(e), u(e), u(e), 1.});
    }
    
    std::vector<double> projNum;
    std::vector<double> counNum;

    result_interpret(s, Resolution, SimBoxL, phi, p0, projNum, p.size());
    count_in_sphere(Radius, SimBoxL, p, p0, counNum);

    double temp0{0}, temp1{0};
    for(auto x : counNum) temp0 += x;
    for(auto x : projNum) temp1 += x;

    std::cout << "Check: " << std::endl;
    std::cout << "-Num_Test_Points = " << TESTPOINTS << "\n";
    std::cout << "-Direct Counting = " << temp0/counNum.size() << "\n";
    std::cout << "-MRA_CS Estimate = " << temp1/projNum.size() << "\n";

    std::cout << "List: " << "\n";
    std::cout << "Num +----[Centre of Sphere]----+ [Counting]  :  [MRA_CS]    :   [delta]" << "\n";
    std::setprecision(6);
    for(int i = 0; i < TESTPOINTS; ++i)
        std::cout << std::setw(3) << i << " | "
                  << std::setw(8) << p0[i].x << std::setw(8) << p0[i].y << std::setw(8) << p0[i].z << " | "
                  << std::setw(9) << counNum[i] << std::setw(14) << projNum[i]
                  << std::setw(15) << projNum[i]-counNum[i] << std::endl;

    
    
}





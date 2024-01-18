#include"mracs.h"

#define NUMTEST  100

int main()
{
    read_parameter();
    auto p = read_in_DM_3vector(DataDirec);
    auto s = sfc(p);
    auto w = wfc(Radius, 0);
    auto c = convol3d(s, w, false);

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, SimBoxL);

    std::vector<Particle> p0;

    for(int i = 0; i < NUMTEST; ++i)
    {   
        p0.push_back(Particle{u(e), u(e), u(e), 1.});
    }
    
    //result_interpret(c, p0, projNum);
    auto projNum = project_value(c,p0,false);
    auto counNum = count_in_sphere(Radius, p, p0);

    const double volume{4./3*M_PI*pow(Radius/SimBoxL*GridLen,3)};
    double temp0{0}, temp1{0}, var{0};;
    for(size_t i = 0; i < p0.size(); ++i) temp0 += counNum[i];
    for(size_t i = 0; i < p0.size(); ++i) temp1 += projNum[i] * volume;
    for(size_t i = 0; i < p0.size(); ++i) var += pow(projNum[i]*volume-counNum[i],2);
    std::cout << "Check: " << std::endl;
    std::cout << "-Num_Test_Points = " << NUMTEST << "\n";
    std::cout << "-Direct Counting = " << temp0/p0.size() << "\n";
    std::cout << "-MRA_CS Estimate = " << temp1/p0.size() << "\n";
    std::cout << "-Var {<delta^2>} = " << var/(p0.size() - 1) << "\n";
    
    
    std::cout << "List: " << "\n";
    std::cout << "Num +----[Centre of Sphere]----+ [Counting]  :  [MRA_CS]    :   [delta]" << "\n";
    std::setprecision(6);
    for(int i = 0; i < NUMTEST; ++i)
        std::cout << std::setw(3) << i << " | "
                  << std::setw(8) << p0[i].x << std::setw(8) << p0[i].y << std::setw(8) << p0[i].z << " | "
                  << std::setw(9) << counNum[i] << std::setw(14) << projNum[i]*volume
                  << std::setw(15) << projNum[i]*volume-counNum[i] << std::endl;

}
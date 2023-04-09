// cylinder check
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_DM_3vector(DataDirec);
    force_kernel_type(4);

    double H{10};
    auto s = sfc(p);
    auto w = wfc(Radius,H);
    auto c = convol3d(s, w);

    auto p0 = generate_random_particle(6,SimBoxL,50);
    auto prj = project_value(c,p0);

    auto cic = count_in_cylinder(Radius,H,p,p0);

    const double norm = M_PI * pow(Radius / SimBoxL * GridLen, 2) * H / SimBoxL * GridLen;
    for(size_t i = 0; i < p0.size(); ++i) prj[i] *= norm;

    double temp0{0}, temp1{0};
    for(size_t i = 0; i < p0.size(); ++i) {
        temp0 += cic[i];
        temp1 += prj[i];
    }

    std::cout << "Check: " << std::endl;
    std::cout << "-Num_Test_Points = " << p0.size() << "\n";
    std::cout << "-Direct Counting = " << temp0/p0.size() << "\n";
    std::cout << "-MRA_CS Estimate = " << temp1/p0.size() << "\n";

    std::cout << "List: " << "\n";
    std::cout << "Num +----[Centre of Sphere]----+ [Counting]  :  [MRA_CS]    :   [delta]" << "\n";
    std::setprecision(6);
    for(int i = 0; i < p0.size(); ++i)
        std::cout << std::setw(3) << i << " | "
                  << std::setw(8) << p0[i].x << std::setw(8) << p0[i].y << std::setw(8) << p0[i].z << " | "
                  << std::setw(9) << cic[i] << std::setw(14) << prj[i]
                  << std::setw(15) << prj[i]-cic[i] << std::endl;

}
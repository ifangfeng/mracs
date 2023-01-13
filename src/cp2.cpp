#include"mracs.h"

#define NUMTEST  50

int main()
{
    read_parameter();
    std::vector<Particle> p;
    p.push_back({567.8, 647.3, 398.6, 1.});
    p.push_back({513.8, 634.3, 878.6, 1.});
    p.push_back({262.8, 347.3, 298.6, 1.});
    const double r0 {5};
    const double r1 {30};
    const double dr {(r1 - r0) / NUMTEST};
    std::vector<double> r(NUMTEST); for(int i = 0; i < NUMTEST; ++i) r[i] = r0 + i * dr;
    std::vector<double> sub0_p1_prj(NUMTEST), sub1_p1_prj(NUMTEST), sub2_p1_prj(NUMTEST), sub3_p1_prj(NUMTEST), sub4_p1_prj(NUMTEST), 
                        sub0_p2_prj(NUMTEST), sub1_p2_prj(NUMTEST), sub2_p2_prj(NUMTEST), sub3_p2_prj(NUMTEST), sub4_p2_prj(NUMTEST), 
                        sub0_p3_prj(NUMTEST), sub1_p3_prj(NUMTEST), sub2_p3_prj(NUMTEST), sub3_p3_prj(NUMTEST), sub4_p3_prj(NUMTEST), 
                        sub0_p1_cic(NUMTEST), sub1_p1_cic(NUMTEST), sub2_p1_cic(NUMTEST), sub3_p1_cic(NUMTEST), sub4_p1_cic(NUMTEST), 
                        sub0_p2_cic(NUMTEST), sub1_p2_cic(NUMTEST), sub2_p2_cic(NUMTEST), sub3_p2_cic(NUMTEST), sub4_p2_cic(NUMTEST), 
                        sub0_p3_cic(NUMTEST), sub1_p3_cic(NUMTEST), sub2_p3_cic(NUMTEST), sub3_p3_cic(NUMTEST), sub4_p3_cic(NUMTEST);
    for(int i = 0; i < NUMTEST; ++i) {
                        sub0_p1_prj[i] = 0, sub1_p1_prj[i] = 0, sub2_p1_prj[i] = 0, sub3_p1_prj[i] = 0, sub4_p1_prj[i] = 0, 
                        sub0_p2_prj[i] = 0, sub1_p2_prj[i] = 0, sub2_p2_prj[i] = 0, sub3_p2_prj[i] = 0, sub4_p2_prj[i] = 0, 
                        sub0_p3_prj[i] = 0, sub1_p3_prj[i] = 0, sub2_p3_prj[i] = 0, sub3_p3_prj[i] = 0, sub4_p3_prj[i] = 0, 
                        sub0_p1_cic[i] = 0, sub1_p1_cic[i] = 0, sub2_p1_cic[i] = 0, sub3_p1_cic[i] = 0, sub4_p1_cic[i] = 0, 
                        sub0_p2_cic[i] = 0, sub1_p2_cic[i] = 0, sub2_p2_cic[i] = 0, sub3_p2_cic[i] = 0, sub4_p2_cic[i] = 0, 
                        sub0_p3_cic[i] = 0, sub1_p3_cic[i] = 0, sub2_p3_cic[i] = 0, sub3_p3_cic[i] = 0, sub4_p3_cic[i] = 0;}

    double* sub0_s = new double[GridNum];
    double* sub1_s = new double[GridNum];
    double* sub2_s = new double[GridNum];
    double* sub3_s = new double[GridNum];
    double* sub4_s = new double[GridNum];
    for(size_t i = 0; i < GridNum; ++i) {sub0_s[i] = 0; sub1_s[i] = 0; sub2_s[i] = 0; sub3_s[i] = 0; sub4_s[i] = 0;}
    const int Nfile{1920};
    for(int i = 0; i < Nfile; ++i)
    {
        std::string ifname = "/data0/tmp/.../snap_130/snap_130." + std::to_string(i);
        std::ifstream ifs(ifname, std::ios_base::binary);
        if(!ifs){
            std::cout << "opening file with error! Abort. File id: " << i << std::endl;
            std::terminate();}
        ifs.seekg(4 + 256 + 4);
        int nbyte;
        int npart;
        ifs.read(as_bytes(nbyte),sizeof(nbyte));
        npart = nbyte / 4 / 3;
        float a[3];
        std::vector<Particle> p0,p1,p2,p3,p4;
        for(int i = 0; i < npart; ++i){
            ifs.read(as_bytes(a),sizeof(a));
            p0.push_back({a[0],a[1],a[2],1.});
        }
        for(int i = 0; i < p0.size(); i+=10) p1.push_back({p0[i].x, p0[i].y, p0[i].z, 1.});
        for(int i = 0; i < p1.size(); i+=10) p2.push_back({p1[i].x, p1[i].y, p1[i].z, 1.});
        for(int i = 0; i < p2.size(); i+=10) p3.push_back({p2[i].x, p2[i].y, p2[i].z, 1.});
        for(int i = 0; i < p3.size(); i+=10) p4.push_back({p3[i].x, p3[i].y, p3[i].z, 1.});

        auto stmp0 = sfc(p0);
        auto stmp1 = sfc(p1);
        auto stmp2 = sfc(p2);
        auto stmp3 = sfc(p3);
        auto stmp4 = sfc(p4);

        for(size_t i = 0; i < GridNum; ++i) {
            sub0_s[i] += stmp0[i];
            sub1_s[i] += stmp1[i];
            sub2_s[i] += stmp2[i];
            sub3_s[i] += stmp3[i];
            sub4_s[i] += stmp4[i];}
        for(int i = 0; i < NUMTEST; ++i){
            auto sub0_cictmp = count_in_sphere(r[i],p0,p);
            auto sub1_cictmp = count_in_sphere(r[i],p1,p);
            auto sub2_cictmp = count_in_sphere(r[i],p2,p);
            auto sub3_cictmp = count_in_sphere(r[i],p3,p);
            auto sub4_cictmp = count_in_sphere(r[i],p4,p);
            sub0_p1_cic[i] += sub0_cictmp[0]; sub0_p2_cic[i] += sub0_cictmp[1]; sub0_p3_cic[i] += sub0_cictmp[2];
            sub1_p1_cic[i] += sub1_cictmp[0]; sub1_p2_cic[i] += sub1_cictmp[1]; sub1_p3_cic[i] += sub1_cictmp[2];
            sub2_p1_cic[i] += sub2_cictmp[0]; sub2_p2_cic[i] += sub2_cictmp[1]; sub2_p3_cic[i] += sub2_cictmp[2];
            sub3_p1_cic[i] += sub3_cictmp[0]; sub3_p2_cic[i] += sub3_cictmp[1]; sub3_p3_cic[i] += sub3_cictmp[2];
            sub4_p1_cic[i] += sub4_cictmp[0]; sub4_p2_cic[i] += sub4_cictmp[1]; sub4_p3_cic[i] += sub4_cictmp[2];
}

        std::cout << "=========>Finshed snap_130." << i << std::endl;
    }

    force_kernel_type(1);
    auto sc0 = sfc_r2c(sub0_s);
    auto sc1 = sfc_r2c(sub1_s);
    auto sc2 = sfc_r2c(sub2_s);
    auto sc3 = sfc_r2c(sub3_s);
    auto sc4 = sfc_r2c(sub4_s);

    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r[i],0);
        double volume = 4./3 * M_PI * pow(r[i],3);
        auto c0 = convol_c2r(sc0,w); auto temp0 = project_value(c0,p);
        auto c1 = convol_c2r(sc1,w); auto temp1 = project_value(c1,p);
        auto c2 = convol_c2r(sc2,w); auto temp2 = project_value(c2,p);
        auto c3 = convol_c2r(sc3,w); auto temp3 = project_value(c3,p);
        auto c4 = convol_c2r(sc4,w); auto temp4 = project_value(c4,p);
        sub0_p1_prj[i] = temp0[0] / volume; sub0_p2_prj[i] = temp0[1] / volume; sub0_p3_prj[i] = temp0[2] / volume;
        sub1_p1_prj[i] = temp0[0] / volume; sub1_p2_prj[i] = temp0[1] / volume; sub1_p3_prj[i] = temp0[2] / volume;
        sub2_p1_prj[i] = temp0[0] / volume; sub2_p2_prj[i] = temp0[1] / volume; sub2_p3_prj[i] = temp0[2] / volume;
        sub3_p1_prj[i] = temp0[0] / volume; sub3_p2_prj[i] = temp0[1] / volume; sub3_p3_prj[i] = temp0[2] / volume;
        sub4_p1_prj[i] = temp0[0] / volume; sub4_p2_prj[i] = temp0[1] / volume; sub4_p3_prj[i] = temp0[2] / volume;
        delete[] w;
        delete[] c0; delete[] c1; delete[] c2; delete[] c3; delete[] c4;
    }

    for(int i = 0; i < NUMTEST; ++i) {
        double volume = 4./3 * M_PI * pow(r[i],3);
        sub0_p1_cic[i] /= volume; sub0_p2_cic[i] /= volume; sub0_p3_cic[i] /= volume;
        sub1_p1_cic[i] /= volume; sub1_p2_cic[i] /= volume; sub1_p3_cic[i] /= volume;
        sub2_p1_cic[i] /= volume; sub2_p2_cic[i] /= volume; sub2_p3_cic[i] /= volume;
        sub3_p1_cic[i] /= volume; sub3_p2_cic[i] /= volume; sub3_p3_cic[i] /= volume;
        sub4_p1_cic[i] /= volume; sub4_p2_cic[i] /= volume; sub4_p3_cic[i] /= volume;
    }

    const double rhobar0 = array_sum(sub0_s,GridNum) / pow(SimBoxL,3);
    const double rhobar1 = array_sum(sub1_s,GridNum) / pow(SimBoxL,3);
    const double rhobar2 = array_sum(sub2_s,GridNum) / pow(SimBoxL,3);
    const double rhobar3 = array_sum(sub3_s,GridNum) / pow(SimBoxL,3);
    const double rhobar4 = array_sum(sub4_s,GridNum) / pow(SimBoxL,3);
    for(auto i : r) std::cout << i / SimBoxL * GridLen << ", " << "======>cic p1 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p1_cic) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p1_cic) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p1_cic) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p1_cic) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p1_cic) std::cout << i / rhobar4 << ", "; std::cout << std::endl << "======>prj p1 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p1_prj) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p1_prj) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p1_prj) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p1_prj) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p1_prj) std::cout << i / rhobar4 << ", "; std::cout << std::endl << "======>cic p2 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p2_cic) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p2_cic) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p2_cic) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p2_cic) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p2_cic) std::cout << i / rhobar4 << ", "; std::cout << std::endl << "======>prj p2 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p2_prj) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p2_prj) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p2_prj) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p2_prj) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p2_prj) std::cout << i / rhobar4 << ", "; std::cout << std::endl << "======>cic p3 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p3_cic) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p3_cic) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p3_cic) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p3_cic) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p3_cic) std::cout << i / rhobar4 << ", "; std::cout << std::endl << "======>prj p3 sub0,1,2,3,4: " << std::endl;
    for(auto i : sub0_p3_prj) std::cout << i / rhobar0 << ", "; std::cout << std::endl;
    for(auto i : sub1_p3_prj) std::cout << i / rhobar1 << ", "; std::cout << std::endl;
    for(auto i : sub2_p3_prj) std::cout << i / rhobar2 << ", "; std::cout << std::endl;
    for(auto i : sub3_p3_prj) std::cout << i / rhobar3 << ", "; std::cout << std::endl;
    for(auto i : sub4_p3_prj) std::cout << i / rhobar4 << ", "; std::cout << std::endl;

}

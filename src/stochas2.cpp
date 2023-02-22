#include"mracs.h"
#include"kdtree.hpp"

#define NUMRAN  1000000
void do_something(std::string cic_or_prj, double* dtm, double* dth, double R);
void write_to_file(std::ofstream& ofs, double* dt, double expect);
void uniform_dt_array(double* dt, double expect);

int main()
{
    read_parameter();
    auto p1 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub005.bin");
    auto p2 = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub2halo.bin");

    std::cout << "dm: "   << p1.size() << std::endl;
    std::cout << "halo: " << p2.size() << std::endl;
    //#################################################

    std::vector<Particle> p0;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int nps = 100;
    const double safeband = 50;

    for(int i = 0; i < nps; ++i)
        for(int j = 0; j < nps; ++j)
            for(int k = 0; k < nps; ++k)
                p0.push_back({safeband + (i + u(e)) * (SimBoxL - 2*safeband) / nps,
                              safeband + (j + u(e)) * (SimBoxL - 2*safeband) / nps,
                              safeband + (k + u(e)) * (SimBoxL - 2*safeband) / nps});

    //#################################################
    double exp_m = p1.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3); 
    double exp_h = p2.size() * 4./3 * M_PI * pow(Radius / SimBoxL,3); 
    std::cout << "exp_m: " << exp_m << std::endl;
    std::cout << "exp_h: " << exp_h << std::endl;

    std::string ofname_prj_dm = "output/prj_R" + std::to_string((int)Radius) + "J" + std::to_string(Resolution) + "_dm.txt";
    std::string ofname_prj_halo  = "output/prj_R" + std::to_string((int)Radius) + "J" + std::to_string(Resolution) + "_halo.txt";
    std::ofstream ofs_prj_dm {ofname_prj_dm};
    std::ofstream ofs_prj_halo {ofname_prj_halo};

    auto w = wfc(Radius, 0);
    auto s1 = sfc(p1); 
    auto s2 = sfc(p2); 
    auto c1 = convol3d(s1,w); 
    auto c2 = convol3d(s2,w);
    auto dtmprj = project_value(c1,p0); delete[] s1; delete[] c1;
    auto dthprj = project_value(c2,p0); delete[] s2; delete[] c2; delete[] w;
    std::cout << "================================prj_R" << Radius << "J" << Resolution << "================: " << std::endl;
    uniform_dt_array(dthprj, exp_h);
    uniform_dt_array(dtmprj, exp_m);
    do_something("prj vs. prj", dtmprj, dthprj, Radius);
    write_to_file(ofs_prj_dm, dtmprj, exp_m);
    write_to_file(ofs_prj_halo, dthprj, exp_h);

    // ==================================================
    // construction of kd-tree
    Kdtree::KdNodeVector nodes;
    double diff;
    clock_t begin, end;
    for (size_t i = 0; i < p1.size(); ++i) {
      std::vector<double> point(3);
      point[0] = p1[i].x;
      point[1] = p1[i].y;
      point[2] = p1[i].z;
      nodes.push_back(Kdtree::KdNode(point));
    }
    begin = clock();
    Kdtree::KdTree tree(&nodes);
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Creation time for " << p1.size() << " p1 dm points:  " << diff << "s" << std::endl;

    // 1.2) range query
    Kdtree::KdNodeVector result;
    auto dtmcic = new double[p0.size()];
    begin = clock();
    for (size_t i = 0; i < p0.size(); ++i) {
      std::vector<double> test_point(3);
      test_point[0] = p0[i].x;
      test_point[1] = p0[i].y;
      test_point[2] = p0[i].z;
      tree.range_nearest_neighbors(test_point, Radius, &result);
      dtmcic[i] = result.size();
    }
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "counting time(tree) for " << p0.size() << " random points:  " << diff << "s" << std::endl;

    // ==================================================
    // construction of kd-tree
    Kdtree::KdNodeVector nodes2;
    for (size_t i = 0; i < p2.size(); ++i) {
      std::vector<double> point2(3);
      point2[0] = p2[i].x;
      point2[1] = p2[i].y;
      point2[2] = p2[i].z;
      nodes2.push_back(Kdtree::KdNode(point2));
    }
    begin = clock();
    Kdtree::KdTree tree2(&nodes2);
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Creation time for " << p2.size() << " p2 dm points:  " << diff << "s" << std::endl;

    // 1.2) range query
    Kdtree::KdNodeVector result2;
    auto dthcic = new double[p0.size()];
    begin = clock();
    for (size_t i = 0; i < p0.size(); ++i) {
      std::vector<double> test_point(3);
      test_point[0] = p0[i].x;
      test_point[1] = p0[i].y;
      test_point[2] = p0[i].z;
      tree2.range_nearest_neighbors(test_point, Radius, &result2);
      dthcic[i] = result2.size();
    }
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "counting time(tree) for " << p0.size() << " random points:  " << diff << "s" << std::endl;

    std::string ofname_cic_dm   = "output/cic_R" + std::to_string((int)Radius) + "_dm.txt"; 
    std::string ofname_cic_halo = "output/cic_R" + std::to_string((int)Radius) + "_halo.txt"; 
    std::ofstream ofs_cic_dm   {ofname_cic_dm};
    std::ofstream ofs_cic_halo {ofname_cic_halo};
    
    std::cout << "++++++++++++++++++++++++++++++++++++cic_R" << Radius << "++++++++++++++++: " << std::endl;
    uniform_dt_array(dthcic, exp_h);
    uniform_dt_array(dtmcic, exp_m);
    do_something("cic vs. mracs", dtmprj, dthcic, Radius);
    do_something("mracs vs. cic", dtmcic, dthprj, Radius);
    do_something("cic vs. cic", dtmcic, dthcic, Radius);
    write_to_file(ofs_cic_dm, dtmcic, exp_m);
    write_to_file(ofs_cic_halo, dthcic, exp_h);
    
}

void write_to_file(std::ofstream& ofs, double* dt, double expect){
    for(size_t i = 0; i < NUMRAN; ++i) ofs << dt[i] << ", ";
}
void uniform_dt_array(double* dt, double expect)
{
    for(size_t i = 0; i < NUMRAN; ++i) dt[i] = dt[i]/expect - 1;
}

void do_something(std::string cic_or_prj, double* dtm, double* dth, double R){

    const int num_bin {15};
    double dtm0 {-1}, dtm1 {2};
    if(R < 20){
        dtm0 = -1;
        dtm1 = 2;
    }
    else if(R < 30){
        dtm0 = -1;
        dtm1 = 1.5;
    }
    else if(R < 50){
        dtm0 = -0.5;
        dtm1 = 1;
    }
    else if(R <= 80){
        dtm0 = -0.4;
        dtm1 = 0.5;
    }
    

    const double ddt {(dtm1 - dtm0) / num_bin};

    std::vector<unsigned> count(num_bin);
    std::vector<double> ave(num_bin), var(num_bin), cbin(num_bin);
    for(int i = 0; i < num_bin; ++i) {ave[i] = 0; var[i] = 0; count[i] = 0;}
    for(int i = 0; i < num_bin; ++i) {cbin[i] = dtm0 + (i + 0.5) * ddt;}

    for(size_t i = 0; i < NUMRAN; ++i){
        int index = floor((dtm[i] - dtm0) / ddt);
        if(index < num_bin && index >= 0){
            ave[index] += dth[i];
            var[index] += pow(dth[i],2);
            ++count[index];
        }
    }

    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]) {
            ave[i] /= count[i];
            var[i] /= count[i];
            var[i] -= pow(ave[i],2);
        }
    }
    double x_tmp{0}, y_tmp{0}, bias1,bias2;
    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]){
            x_tmp += cbin[i];
            y_tmp += ave[i];
        }
    }
    bias1 = y_tmp/x_tmp;
    x_tmp = 0;
    y_tmp = 0;
    for(size_t i = 0; i < num_bin; ++i){
        if(count[i]){
            x_tmp += cbin[i] * count[i];
            y_tmp += ave[i] * count[i];
        }
    }
    bias2 = y_tmp/x_tmp;
    std::cout << "bias1= " << bias1 << std::endl;
    std::cout << "bias2= " << bias2 << std::endl;
    std::cout << "centre of bin: " << std::endl; for(auto i : cbin) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "count: " << std::endl; for(auto i : count) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h: " << std::endl; for(auto i : ave) std::cout << i << ", "; std::cout << std::endl;
    std::cout << "delta_h deviation: " << std::endl; for(auto i : var) std::cout << sqrt(i) << ", "; std::cout << std::endl;

    std::cout << "=================>finished one sub map." << std::endl;
}


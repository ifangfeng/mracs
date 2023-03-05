// noise cancelling test
#include"mracs.h"
#include"kdtree.hpp"
using namespace std;


int main()
{
    read_parameter();
    std::string fname{"dm_sub5e-6"};
    auto p = read_in_DM_3vector("/data0/MDPL2/dm_sub/" + fname + ".bin");
    // std::string fname{"halo_position"};
    // auto p1 = read_in_Halo_4vector("/data0/MDPL2/" + fname + ".bin");
    // std::vector<Particle> p;
    // for(auto x : p1) if(x.weight > 2e12) p.push_back({x.x,x.y,x.z,x.weight});
    // std::cout << "number of halo: " << p.size() << std::endl;
    // std::vector<Particle>().swap(p1);

    auto p0 = generate_random_particle(215,SimBoxL,50);

    clock_t begin, end;
    double diff;
   
     // 1.1) construction of kd-tree
    Kdtree::KdNodeVector nodes;

    for (size_t i = 0; i < p.size(); ++i) {
      std::vector<double> point(3);
      point[0] = p[i].x;
      point[1] = p[i].y;
      point[2] = p[i].z;
      nodes.push_back(Kdtree::KdNode(point));
    }
    begin = clock();
    Kdtree::KdTree tree(&nodes);
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Creation time for " << p.size() << " galaxy points:  " << diff << "s" << endl;

    // 1.2) range query
    Kdtree::KdNodeVector result;
    std::vector<int64_t> count(p0.size());
    begin = clock();
    for (size_t i = 0; i < p0.size(); ++i) {
      std::vector<double> test_point(3);
      test_point[0] = p0[i].x;
      test_point[1] = p0[i].y;
      test_point[2] = p0[i].z;
      tree.range_nearest_neighbors(test_point, Radius, &result);
      count[i] = result.size();
    }
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    cout << "counting time(tree) for " << p0.size() << " random points:  " << diff << "s" << endl;

    double n = 4./3 * M_PI * pow(Radius/SimBoxL,3) * p.size();
    cic_pdf(count, 0, 5, n, fname);

    
}



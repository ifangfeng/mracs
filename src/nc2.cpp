// noise cancelling test
#include"mracs.h"
#include"kdtree.hpp"
using namespace std;

int main()
{
    read_parameter();
    std::vector<std::string> namestr = {"05", "005", "5e-4", "5e-5", "2halo", "5e-6"};
    auto p = read_in_DM_3vector(DataDirec);
    auto p0 = generate_random_particle(10,SimBoxL,50);

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
    
    return 0;
}
#include"mracs.h"
#include"kdtree.hpp"


using namespace std;

//
// helper function for printing points
//
void print_nodes(const Kdtree::KdNodeVector &nodes) {
  size_t i,j;
  for (i = 0; i < nodes.size(); ++i) {
    if (i > 0)
      cout << " ";
    cout << "(";
    for (j = 0; j < nodes[i].point.size(); j++) {
      if (j > 0)
        cout << ",";
      cout << nodes[i].point[j];
    }
    cout << ")";
  }
  cout << endl;
}

//
// main program demonstrating typical use cases
//
int main(int argc, char** argv) {

  read_parameter();
  auto p = read_in_DM_3vector(DataDirec);

  double diff;
  clock_t begin, end;

  // 0.0 generating random sampling points
  std::vector<Particle> p0;
  std::default_random_engine e;
  std::uniform_real_distribution<double> u(0,1);
  const int nps = 5;
  const double safeband = 50;

  begin = clock();
  for(int i = 0; i < nps; ++i)
      for(int j = 0; j < nps; ++j)
          for(int k = 0; k < nps; ++k)
              p0.push_back({safeband + (i + u(e)) * (SimBoxL - 2*safeband) / nps,
                            safeband + (j + u(e)) * (SimBoxL - 2*safeband) / nps,
                            safeband + (k + u(e)) * (SimBoxL - 2*safeband) / nps});
  end = clock();
  diff = double(end - begin) / CLOCKS_PER_SEC;
  cout << "generating time for " << p0.size() << " random points:  " << diff << "s" << endl;


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

  // 2.1 count in sphere
  auto cic = count_in_sphere(Radius, p, p0);

  // 3.1 check
  cout << "tree count: " << endl;
  for(int i = 0; i < p0.size(); ++i)
    cout << count[i] << ", ";
  cout << endl << "count in sphere: " << endl;
  for(int i = 0; i < p0.size(); ++i)
    cout << cic[i] << ", ";
  cout << endl;
  return 0;
}

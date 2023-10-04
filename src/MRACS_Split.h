#include"MRACS_Main.h"

// which vector should trial been push back
template<class T> int classify_index(std::vector<T>& node, T trial){
    int index{0};
    for(int i = 0; i < node.size(); ++i){
        if(trial > node[i]) ++index;
        else break;
    }
    return index;
}

Point eigen_element(double xx, double xy, double xz, double yy, double yz, double zz);
std::vector<Point> eigenvalue_of_tidal_tensor(double** cxx,std::vector<Particle>& p);
double** tidal_tensor(fftw_complex* sc, double* w);
double* tensor_element(fftw_complex* sc, uint dim_i, uint dim_j);
std::vector<std::vector<Particle>*> halo_mass_split(std::vector<Particle>& hl, int nbin);
void print_min_max_and_size(std::vector<Particle>& hl);
std::vector<int> envi_with_Mcut(std::string ifn, double Mcut, std::vector<Particle>& hl);
void print_min_max_and_size_double(std::vector<double>& vec);
std::vector<std::vector<Particle>*> halo_envi_match_and_split(std::vector<int>& envi, std::vector<Particle>& hl);
std::vector<int> envi_vector_readin(std::string ifn, size_t hlsize);
std::vector<std::vector<Particle>*> mass_classify_and_push_back(std::vector<double>& node, std::vector<Particle>& hl, int nbin);
std::vector<std::vector<Particle>*> halo_envi_mass_multi_split(std::vector<int>& envi, std::vector<Particle>& hl, int nbin);
std::vector<std::vector<Particle>*> halo_envi_mass_concatenate_split(std::vector<int>& envi, std::vector<Particle>& hl, int nbin);
std::vector<double> nodes_of_proto_sort(std::vector<double>& vec, int nbin);
std::vector<size_t> limited_sort(std::vector<double>& vec);
size_t minimum_index(std::vector<double>& v);
size_t maximum_index(std::vector<double>& v);
double** tidal_tensor(fftw_complex* sc, double* w);
int eigen_classify(double xx, double xy, double xz, double yy, double yz, double zz, double lambda_th);
std::vector<int> web_classify(double** cxx, std::vector<Particle>& p, double lambda_th);
std::vector<int> web_classify_to_grid(double** cxx, double lambda_th);
std::vector<int> environment(std::vector<Particle>& dm, double Rs, std::vector<Particle>& p0);
double gaussian_radius_from_mass(double m_smooth);
fftw_complex* hermitian_product(fftw_complex* sc1, fftw_complex* sc2);
std::vector<double> optimal_weight_solver(std::vector<double> cov, int n, bool PRINT);
fftw_complex* optimal_reconstruct(fftw_complex* sc_dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT);
std::vector<double> optimal_solution(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT);
std::vector<double> optimal_solution_lean(fftw_complex* sc_dm, std::vector<std::vector<Particle>*> vpts, double* wpk, bool PRINT);
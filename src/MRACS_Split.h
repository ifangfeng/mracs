#include"MRACS_Main.h"

std::vector<std::vector<Particle>*> halo_mass_split(std::vector<Particle>& hl, int nbin);
void print_min_max_and_size(std::vector<Particle>& hl);
void print_min_max_and_size_double(std::vector<double>& vec);
int classify_index(std::vector<double>& node, double trial);
std::vector<int> envi_vector_readin(std::string ifn, size_t hlsize);
std::vector<std::vector<Particle>*> mass_classify_and_push_back(std::vector<double>& node, std::vector<Particle>& hl, int nbin);
std::vector<std::vector<Particle>*> halo_envi_mass_multi_split(std::string ifn, std::vector<Particle>& hl, int nbin);
std::vector<std::vector<Particle>*> halo_envi_mass_concatenate_split(std::string ifn, std::vector<Particle>& hl, int nbin);
std::vector<double> nodes_of_proto_sort(std::vector<double>& vec, int nbin);
std::vector<size_t> limited_sort(std::vector<double>& vec);
size_t minimum_index(std::vector<double>& v);
size_t maximum_index(std::vector<double>& v);
double** tidal_tensor(fftw_complex* sc, double* w);
int eigen_classify(double xx, double xy, double xz, double yy, double yz, double zz);
std::vector<int> web_classify(double** cxx, std::vector<Particle>& p);
std::vector<int> web_classify_to_grid(double** cxx);
std::vector<int> environment(std::vector<Particle>& dm, double Rs, std::vector<Particle>& p0);
double gaussian_radius_from_mass(double m_smooth);
fftw_complex* hermitian_product(fftw_complex* sc1, fftw_complex* sc2);
std::vector<double> optimal_weight_solver(std::vector<double> cov, int n, bool PRINT);
std::vector<std::vector<Particle>*> halo_envi_match_and_split(std::string ifn, std::vector<Particle>& hl);
fftw_complex* optimal_reconstruct(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT);
std::vector<double> optimal_solution(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT);

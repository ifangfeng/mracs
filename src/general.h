#include"csmain.h"

std::vector<double> log_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT);
std::vector<double> linear_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT);
std::vector<Particle> default_random_particle(double boxsize, size_t n);
std::vector<Particle> generate_random_particle(int x, double L, double w);
void fill_index_set(const double R, std::vector<Index>& inner_index, std::vector<Index>& cross_index);
double* count_in_sphere(const double R, std::vector<Particle>& p, std::vector<Particle>& p0);
double* count_in_cylinder(double R, double H, std::vector<Particle>& p, std::vector<Particle>& p0);

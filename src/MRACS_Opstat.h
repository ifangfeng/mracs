#include"MRACS_Main.h"


void pdf(std::vector<Particle>& p0, double* c, double nf, double rhomin, double rhomax, int nbin, std::string ofname);
void cp_dispersion(std::vector<int64_t>& c, double* n, double rhomin, double rhomax, double cicexpect, std::string ofname);
void cic_pdf(std::vector<int64_t>& c, double rhomin, double rhomax, double cicexpect, std::string ofname);

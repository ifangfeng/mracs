// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl);

int main(){
    read_parameter();
    
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_3vector("/data0/MDPL2/halo_Mcut2e12.bin");

    auto ccc = fourier_mode_correlation_1rlz(dm,hl);

    std::ofstream ofs {"output/data.txt"}; for(auto x : ccc) ofs << x << " ";
    std::ofstream ofs1 {"output/datax.txt"}; for(int i = 0; i < GridLen/2 + 1; ++i) ofs1 << i * TWOPI / SimBoxL << " ";
}



std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl)
{
    auto dm_sc = sfc_r2c(sfc(dm),true);
    auto hl_sc = sfc_r2c(sfc(hl),true);

    double mh[GridLen](0);
    double mm[GridLen](0);
    double hh[GridLen](0);
    
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < (GridLen/2 + 1); ++k){
                auto ii = i < GridLen/2 + 1 ? i : GridLen - i;
                auto jj = j < GridLen/2 + 1 ? j : GridLen - j;
                int ll = sqrt(ii * ii + jj * jj + k * k);
                auto l = i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k;

                mh[ll] += dm_sc[l][0] * hl_sc[l][0] + dm_sc[l][1] * hl_sc[l][1];
                mm[ll] += dm_sc[l][0] * dm_sc[l][0] + dm_sc[l][1] * dm_sc[l][1];
                hh[ll] += hl_sc[l][0] * hl_sc[l][0] + hl_sc[l][1] * hl_sc[l][1];
            }
    std::vector<double> cross(GridLen/2 + 1);
    for(int i = 0; i < GridLen/2 + 1; ++i) cross[i] = mh[i]/sqrt(mm[i] * hh[i]);
    fftw_free(dm_sc);
    fftw_free(hl_sc);

    return cross;
}
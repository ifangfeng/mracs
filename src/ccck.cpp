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

std::vector<std::vector<Particle>> density_split(std::vector<Particle>& hl, double R, int T, int Nbin)
{

}

std::vector<std::vector<Particle>> halo_mass_split(std::vector<Particle>& hl, int nbin)
{
    size_t max_id{0}, min_id{0};

    for(size_t i = 0; i < hl.size(); ++i){
        if (hl[i].weight > hl[max_id].weight) 
            max_id = i;
        else if(hl[i].weight < hl[min_id].weight)
            min_id = i;
    }
    const size_t MAXID{maximum_index(hl)}, MINID{min_id};
    const double DELTA_MASS{hl[MAXID].weight - hl[MINID].weight};
    size_t dx = hl.size() / nbin;
    const int INITIAL{10}, REFINE{1000};
    const size_t LINSIZE{hl.size() < 1e10 ? INITIAL * hl.size() : hl.size()};
    int* count = new int[LINSIZE](0);

    for(size_t i = 0; i < hl.size(); ++i){
        size_t idx = (hl[i].weight - hl[MINID].weight)/DELTA_MASS * LINSIZE;
        ++count[idx];
    }
    double node[nbin-1](0);
    size_t tmp{0};
    size_t sum{0};
    for(size_t i = 0; i < LINSIZE; ++i){
        sum += count[i];
        if(sum >= hl.size()/nbin){
            tmp = i;
            break;
        }
    }
    std::vector<double> nodePart;
    for(size_t i = 0; i < hl.size(); ++i){
        size_t idx = (hl[i].weight - hl[MINID].weight)/DELTA_MASS * LINSIZE;
        if(idx = tmp) nodePart.push_back(hl[i].weight);
    }
}

std::vector<size_t> proto_sort(std::vector<double>& vec, int nbin)
{
    auto max_id = maximum_index(vec);
    auto min_id = minimum_index(vec);

    double MAX{vec[max_id]}, MIN{vec[min_id]};
    double DELTA{MAX - MIN};
    const int REFINE{100};
    const size_t LINSIZE{vec.size() < 1e9 ? vec.size() * REFINE : vec.size()};

    auto count = new int[LINSIZE](0);
    for(size_t i = 0; i < vec.size(); ++i){
        size_t idx = (vec[i] - MIN) / DELTA * LINSIZE;
        ++count[idx];
    }

    double sum{0};
    size_t node_id{0}, finer_id{0};
    const size_t dx {vec.size()/nbin};

    for(size_t i = 0; i < LINSIZE; ++i){
        sum += count[i];
        if(sum >= dx) 
        {   
            node_id = i;
            break;
        }
    }
    finer_id =  dx - (sum - count[node_id]);
}

// -------------------------------------------------------
// direct sort of double vector, return as ascend index
// notice the original vector will be modified
// -------------------------------------------------------
std::vector<size_t> limited_sort(std::vector<double>& vec)
{
    std::vector<size_t> sortedID;
    size_t max_id = maximum_index(vec);
    const double MAX{vec[max_id]};
    for(size_t i = 0; i < vec.size() - 1; ++i){
        size_t id = minimum_index(vec);
        vec[id] = MAX;
        sortedID.push_back(id);
    }
    sortedID.push_back(max_id);

    return sortedID;
}

size_t minimum_index(std::vector<double>& v)
{
    size_t idx{0};
    for(size_t i = 0; i < v.size(); ++i){
        if(v[i] < v[idx]) idx = i;
    }
    return idx;
}

size_t maximum_index(std::vector<double>& v)
{
    size_t idx{0};
    for(size_t i = 0; i < v.size(); ++i){
        if(v[i] > v[idx]) idx = i;
    }
    return idx;
}
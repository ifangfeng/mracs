// *******************************************************
// halo environment split and cross-correlation function
// *******************************************************
#include"mracs.h"

std::vector<double> fourier_mode_correlation_1rlz(std::vector<Particle>& dm, std::vector<Particle>& hl);
void print_min_max_and_size(std::vector<Particle>& hl);
std::vector<std::vector<Particle>> halo_mass_split(std::vector<Particle>& hl, int nbin);
int classify_index(std::vector<double>& node, double trial);
std::vector<multipoint> proto_sort(std::vector<double>& vec, int nbin);
std::vector<size_t> limited_sort(std::vector<double> vec);
size_t minimum_index(std::vector<double>& v);
size_t maximum_index(std::vector<double>& v);

int main(){
    read_parameter();
    
    //auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    auto hl = read_in_Halo_4vector("/data0/MDPL2/halo_Mcut2e12.bin");
    auto cata = halo_mass_split(hl,5);

    for(auto x : cata) print_min_max_and_size(x);
    // auto ccc = fourier_mode_correlation_1rlz(dm,hl);
    // std::ofstream ofs {"output/data.txt"}; for(auto x : ccc) ofs << x << " ";
    // std::ofstream ofs1 {"output/datax.txt"}; for(int i = 0; i < GridLen/2 + 1; ++i) ofs1 << i * TWOPI / SimBoxL << " ";
}

void print_min_max_and_size(std::vector<Particle>& hl){
    double min{hl[0].weight},max{hl[0].weight};
    for(auto x : hl){
        if(x.weight > max) max = x.weight;
        else if(x.weight < min) min = x.weight;
    }
    std::cout << "size: " << hl.size() << ", min: " << min << ", max: " << max << std::endl; 
}



//std::vector<std::vector<Particle>> density_split(std::vector<Particle>& hl, double R, int T, int Nbin)
//{
//
//}

std::vector<std::vector<Particle>> halo_mass_split(std::vector<Particle>& hl, int nbin)
{
    std::vector<double> vecmass(hl.size());
    #pragma omp parallel for
    for(size_t i = 0; i < hl.size(); ++i) vecmass[i] = hl[i].weight;
    auto node = proto_sort(vecmass, nbin);
    for(auto x : node) std::cout << "node: " << x << std::endl;

    int idx{0};
    std::vector<std::vector<Particle>> cata(nbin);
    for(auto x : hl){
        cata[classify_index(node,x.weight)].push_back(x);
    }
    for(auto x : cata) if(x.size() - hl.size()/nbin > nbin) {
        std::cout << "[func: SPLIT] !Warning, some nodes include multiple identical values\n";
        break;
    }

    return cata;
}

// which vector should trial been push back
int classify_index(std::vector<double>& node, double trial){
    int index{0};
    for(int i = 0; i < node.size(); ++i){
        if(trial > node[i]) ++index;
        else break;
    }
    return index;
}

struct multipoint
{
    double x;
    size_t n;
};

// **************************************************************************
// return the nbin fraction node points of a double vector in ascending order
// **************************************************************************
std::vector<multipoint> proto_sort(std::vector<double>& vec, int nbin)
{
    std::vector<multipoint> node;
    auto max_id = maximum_index(vec);
    auto min_id = minimum_index(vec);

    double MAX{vec[max_id]}, MIN{vec[min_id]};
    double DELTA{MAX - MIN};
    const int REFINE{100};
    const size_t LINSIZE{vec.size() < 1e9 ? vec.size() * REFINE : vec.size()};
    const size_t VECLEN{vec.size()/nbin};

    auto count = new int[LINSIZE + 1](0);
    for(size_t i = 0; i < vec.size(); ++i){
        size_t idx = (vec[i] - MIN) / DELTA * LINSIZE;
        ++count[idx];
    }
    
    size_t node_id{0};
    for(int n = 0; n < nbin - 1; ++n){
        size_t sum{0}, finer{0};
        for(size_t i = node_id; i < LINSIZE + 1; ++i){
            sum += count[i];
            if(sum >= VECLEN) 
            {   
                node_id = i;
                break;
            }
        }
        if(count[node_id] > 1){
            finer =  VECLEN - (sum - count[node_id]);
            std::vector<double> nodex;
            for(size_t i = 0; i < vec.size(); ++i){
                size_t idx = (vec[i] - MIN) / DELTA * LINSIZE;
                if(idx == node_id) nodex.push_back(vec[i]);
            }
            auto index = limited_sort(nodex);
            node.push_back(nodex[index[finer-1]]);
            std::vector<double>().swap(nodex);
            std::vector<size_t>().swap(index);
            count[node_id] -= finer;
        }
        else 
            node.push_back(node_id * DELTA / LINSIZE + MIN);
        
    }
    delete count;

    return node;
}

// -------------------------------------------------------
// direct sort of double vector, return as ascending index
// -------------------------------------------------------
std::vector<size_t> limited_sort(std::vector<double> vec)
{
    std::cout << "[func: limited_sort] node size: " << vec.size() << std::endl;
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
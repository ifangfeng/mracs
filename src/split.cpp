// =======================================================
// This cpp source file is part of MRACS project
// halo mass split and environmental classify 
// of cosmic web and related 
// =======================================================

#include"split.h"


// ************************************************************************************
// nbin is the number of catalogue splited, return as vector of catalog
// ************************************************************************************
std::vector<std::vector<Particle>> halo_mass_split(std::vector<Particle>& hl, int nbin)
{
    std::vector<double> vecmass(hl.size());
    #pragma omp parallel for
    for(size_t i = 0; i < hl.size(); ++i) vecmass[i] = hl[i].weight;
    auto node = proto_sort(vecmass, nbin);

    std::vector<std::vector<Particle>> cata(nbin);
    for(auto x : hl){
        cata[classify_index(node,x.weight)].push_back(x);
    }
    for(auto x : cata) if(x.size() - hl.size()/nbin > nbin) {
        std::cout << "[func: SPLIT] !Warning, some nodes include multiple identical items\n";
        break;
    }
    for(auto x : cata) print_min_max_and_size(x);

    return cata;
}

void print_min_max_and_size(std::vector<Particle>& hl){
    double min{hl[0].weight},max{hl[0].weight};
    for(auto x : hl){
        if(x.weight > max) max = x.weight;
        else if(x.weight < min) min = x.weight;
    }
    std::cout << "size: " << hl.size() << ", min: " << min << ", max: " << max << std::endl; 
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

// ************************************************************************************
// return the nbin fraction node points of a double vector in ascending order
// ************************************************************************************
std::vector<double> proto_sort(std::vector<double>& vec, int nbin)
{
    std::vector<double> node;
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

// ************************************************************************************
// direct sort of double vector, return as ascending index
// ************************************************************************************
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






// ************************************************************************************
// given a Cloud-in-Cell interpreted 3d density array in fourier space, and the Gaussian
// kernel also in fourier space, calculate the soothed density fileds and the Hessian 
// of the gravitation potential fileds (df, xx, xy, xz, yy, yz, zz).return as pointer of 
// <double> pointer cxx, with cxx[0] the address of df, cxx[1] address of xx and so on.
// ************************************************************************************
double** tidal_tensor(fftw_complex* sc, double* w)
{
    auto sc0 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    double** cxx = new double*[7];
    for(int i = 0; i < 7; ++i){
        cxx[i] = new double[GridVol];
    }
    fftw_plan pldf = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc0, cxx[0], FFTW_MEASURE);
    fftw_plan plxx = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[1], FFTW_MEASURE);
    fftw_plan plxy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[2], FFTW_MEASURE);
    fftw_plan plxz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[3], FFTW_MEASURE);
    fftw_plan plyy = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[4], FFTW_MEASURE);
    fftw_plan plyz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[5], FFTW_MEASURE);
    fftw_plan plzz = fftw_plan_dft_c2r_3d(GridLen, GridLen, GridLen, sc1, cxx[6], FFTW_MEASURE);

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    {
        sc0[i][0] = w[i] * sc[i][0];
        sc0[i][1] = w[i] * sc[i][1]; 
    }
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * i / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * i / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxx);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * j / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * i * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * i * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plxz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * j / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * j * j / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plyy);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * j * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * j * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plyz);
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 + 1; ++k){
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][0] * k * k / (i*i + j*j + k*k);
                sc1[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] = 
                sc0[i * GridLen * (GridLen/2 + 1) + j * (GridLen/2 + 1) + k][1] * k * k / (i*i + j*j + k*k); 
            }
    sc1[0][0] = 0;
    sc1[0][1] = 0;
    fftw_execute(plzz);
    fftw_execute(pldf);
    
    fftw_destroy_plan(pldf);
    fftw_destroy_plan(plxx);
    fftw_destroy_plan(plxy);
    fftw_destroy_plan(plxz);
    fftw_destroy_plan(plyy);
    fftw_destroy_plan(plyz);
    fftw_destroy_plan(plzz);
    fftw_free(sc0);
    fftw_free(sc1);

    return cxx;
}




// ************************************************************************************
// solving a 3x3 real symmetric matrix and return the number of eigenvalue 
// above a threshold Lambda_th = 0
// ************************************************************************************
int eigen_classify(double xx, double xy, double xz, double yy, double yz, double zz)
{
    int n = 0;
    double lambda_th = 0;
    double t = xx + yy + zz;
    double mid1 = 2 * xx - yy - zz;
    double mid2 = 2 * yy - xx - zz;
    double mid3 = 2 * zz - xx - yy;
    double a = xx * xx + yy * yy + zz * zz - xx * yy - xx * zz - yy * zz + 3 * (xy * xy + xz * xz + yz * yz);
    double b = 9 * (mid1 * yz * yz + mid2 * xz * xz + mid3 * xy * xy) - 54 * xy * xz * yz - mid1 * mid2 * mid3;
    double phi = M_PI / 6;
    if(b > 0)
        phi = atan(sqrt(4 * a * a * a - b * b) / b) / 3.;
    else if(b < 0)
        phi = (atan(sqrt(4 * a * a * a - b * b) / b) + M_PI) / 3.;
    double lambda[3];
    lambda[0] = (t - 2 * sqrt(a) * cos(phi)) / 3;
    lambda[1] = (t - 2 * sqrt(a) * cos(phi + 2 * M_PI / 3)) / 3;
    lambda[2] = (t - 2 * sqrt(a) * cos(phi - 2 * M_PI / 3)) / 3;
    for(int i = 0; i < 3; ++i)
        if(lambda[i] > lambda_th)
            ++n;
    return n;
}


// ************************************************************************************
// from Hessian to web structure: 0-Voids, 1-sheets, 2-filaments, 3-knots
// ************************************************************************************
std::vector<int> web_classify(double** cxx, std::vector<Particle>& p)
{
    std::vector<int> s(p.size());
    #pragma omp parallel for
    for(int i = 0; i < p.size(); ++i)
    {
        int xs,ys,zs;   // BSpline have support [0,n+1],not centre in origin
        int64_t x,y,z,l;   
        xs = p[i].x / SimBoxL * GridLen + 0.5;
        ys = p[i].y / SimBoxL * GridLen + 0.5;
        zs = p[i].z / SimBoxL * GridLen + 0.5;
        x = (xs - 1) & (GridLen - 1);   // shift -1 for CIC, which is corresponding to BSpline n=1
        y = (ys - 1) & (GridLen - 1);
        z = (zs - 1) & (GridLen - 1);
        l = x * GridLen * GridLen + y * GridLen + z;
        s[i] = eigen_classify(cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], cxx[6][l]);
    }
    return s;
}


// ************************************************************************************
// from Hessian to web structure on grid: 0-Voids, 1-sheets, 2-filaments, 3-knots
// ************************************************************************************
std::vector<int> web_classify_to_grid(double** cxx)
{
    std::vector<int> s(GridVol);
    #pragma omp parallel for
    for(int64_t i = 0; i < GridLen; ++i)
        for(int64_t j = 0; j < GridLen; ++j)
            for(int64_t k = 0; k < GridLen; ++k)
            {
                int64_t x,y,z,l;   // BSpline have support [0,n+1],not centre in origin
                x = (i - 1) & (GridLen - 1);   // shift -1 for CIC, which is corresponding to BSpline n=1
                y = (j - 1) & (GridLen - 1);
                z = (k - 1) & (GridLen - 1);
                l = x * GridLen * GridLen + y * GridLen + z;
                s[i * GridLen * GridLen + j * GridLen + k] = eigen_classify(cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], cxx[6][l]);
            }
    return s;
}


// ************************************************************************************
// given dark matter fields dm and Gaussian smoothing radius Rs, return the web classify result at point p0
// ************************************************************************************
std::vector<int> environment(std::vector<Particle>& dm, double Rs, std::vector<Particle>& p0)
{
    force_base_type(0,1);
    force_kernel_type(2);
    auto begin = std::chrono::steady_clock::now();

    auto sc = sfc_r2c(sfc(dm),true);
    auto w = wft(2, 0);
    auto cxx = tidal_tensor(sc, w);
    auto env = web_classify(cxx,p0);
    
    fftw_free(sc);
    delete[] w;

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time difference EnvironmentClassify  = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
    << "[ms]" << std::endl;

    return env;
}


// ************************************************************************************
// enviriamental parameter array locate in gride point, defined as the number of positive eigenvalue
// of tidle tensor of matter density fileds, which is obtand by Cloud-in-Cell interpolation of particles
// to grid point and then smoothed by a Gaussian kernel with radius R. The main process is working on 
// fourier space so we can take advantage of FFT, for detials see Hahn O., Porciani C., Carollo C. M., Dekel A., 2007, MNRAS, 375, 489
// https://ui.adsabs.harvard.edu/abs/2007MNRAS.375..489H
// ************************************************************************************
double gaussian_radius_from_mass(double m_smooth) 
{
    return 1. / sqrt(TWOPI) * pow(m_smooth, 1./3); // rho_bar is needed: [m_smooth / rho_ar]
}

// ************************************************************************************
// p[i] = s1[i] * Hermitian[s2[i]]
// ************************************************************************************
fftw_complex* hermitian_product(fftw_complex* sc1, fftw_complex* sc2)
{
    auto c = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));

    #pragma omp parallel for
    for(size_t i = 0; i < GridLen * GridLen * (GridLen/2 + 1); ++i)
    {
        c[i][0] = sc1[i][0] * sc2[i][0] + sc1[i][1] * sc2[i][1];
        c[i][1] = sc1[i][1] * sc2[i][0] - sc1[i][0] * sc2[i][1];
    }

    return c;
}
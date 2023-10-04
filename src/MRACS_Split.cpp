// =======================================================
// This cpp source file is part of MRACS project
// halo mass split and environmental classify 
// of cosmic web and related 
// =======================================================

#include"MRACS_Split.h"
#include"MRACS_Corr.h"


//============================
double* tensor_element(fftw_complex* sc, uint dim_i, uint dim_j)
{
    if(dim_i > 2 || dim_j > 2){std::cout << "[func: tensor()]: input error!\n"; std::terminate();}

    auto Tk = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    auto T  = new double[GridVol];
    
    #pragma omp parallel for
    for(size_t i = 0; i < GridLen; ++i)
        for(size_t j = 0; j < GridLen; ++j)
            for(size_t k = 0; k < GridLen/2 +1; ++k){
                size_t a[3]{i,j,k};
                double kmodsq = i * i + j * j + k * k;
                Tk[i * GridLen * (GridLen/2 +1) + j * (GridLen/2 + 1) + k][0] =  a[dim_i] * a[dim_j]/kmodsq * 
                sc[i * GridLen * (GridLen/2 +1) + j * (GridLen/2 + 1) + k][0];
                Tk[i * GridLen * (GridLen/2 +1) + j * (GridLen/2 + 1) + k][1] =  a[dim_i] * a[dim_j]/kmodsq * 
                sc[i * GridLen * (GridLen/2 +1) + j * (GridLen/2 + 1) + k][1];
            }Tk[0][0] = 0; Tk[0][1] = 0;

    auto plan = fftw_plan_dft_c2r_3d(GridLen,GridLen,GridLen,Tk,T,FFTW_MEASURE);
    fftw_execute(plan);

    #pragma omp parallel for
    for(size_t i = 0; i < GridVol; ++i) T[i] /= GridVol;

    fftw_destroy_plan(plan);
    fftw_free(Tk);

    return T;
}


// solving normalized Poisson equation in Fourier sapce
double** tidal_tensor(fftw_complex* sc, double* w){
    auto sc1 = fftw_alloc_complex(GridLen * GridLen * (GridLen/2 + 1));
    #pragma omp parallel for
    for(size_t i = 1; i < GridLen * GridLen * (GridLen/2 + 1); ++i){
        sc1[i][0] = sc[i][0] * w[i] / (sc[0][0] * w[0]);
        sc1[i][1] = sc[i][1] * w[i] / (sc[0][0] * w[0]);
    }sc1[0][0] = 0;
    const uint dim{3};

    auto T = new double*[6];
    int n{0};
    for(int i = 0; i < dim; ++i)
        for(int j = i; j < dim; ++j){
            T[n] = tensor_element(sc,i,j);
            ++n;
        }
    return T;
}

// ************************************************************************************
// solving a 3x3 real symmetric matrix and return the number of eigenvalue 
// above a threshold Lambda_th = 0
// ************************************************************************************
int eigen_classify(double xx, double xy, double xz, double yy, double yz, double zz, double lambda_th)
{
    int n = 0;
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

std::vector<Point> eigenvalue_of_tidal_tensor(double** cxx,std::vector<Particle>& p){
    std::vector<Point> vec_eigen(p.size());
    #pragma omp parallel for
    for(size_t i = 0; i < p.size(); ++i){
        int xs,ys,zs;   // BSpline have support [0,n+1],not centre in origin
        int64_t x,y,z,l;   
        xs = p[i].x / SimBoxL * GridLen + 0.5;
        ys = p[i].y / SimBoxL * GridLen + 0.5;
        zs = p[i].z / SimBoxL * GridLen + 0.5;
        x = (xs) & (GridLen - 1);   // shift -1 for CIC, which is corresponding to BSpline n=1
        y = (ys) & (GridLen - 1);
        z = (zs) & (GridLen - 1);
        l = x * GridLen * GridLen + y * GridLen + z;
        vec_eigen[i] = eigen_element(cxx[0][l], cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l]);
    }
    return vec_eigen;
}

Point eigen_element(double xx, double xy, double xz, double yy, double yz, double zz){
    int n = 0;
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

    Point lambda;
    lambda.x = (t - 2 * sqrt(a) * cos(phi)) / 3;
    lambda.y = (t - 2 * sqrt(a) * cos(phi + 2 * M_PI / 3)) / 3;
    lambda.z = (t - 2 * sqrt(a) * cos(phi - 2 * M_PI / 3)) / 3;

    return lambda;
}



// ************************************************************************************
// from Hessian to web structure: 0-Voids, 1-sheets, 2-filaments, 3-knots
// ************************************************************************************
std::vector<int> web_classify(double** cxx, std::vector<Particle>& p, double lambda_th)
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
        s[i] = eigen_classify(cxx[0][l], cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], lambda_th);
    }
    return s;
}


// ************************************************************************************
// from Hessian to web structure on grid: 0-Voids, 1-sheets, 2-filaments, 3-knots
// ************************************************************************************
std::vector<int> web_classify_to_grid(double** cxx, double lambda_th)
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
                s[i * GridLen * GridLen + j * GridLen + k] = eigen_classify(cxx[0][l], cxx[1][l], cxx[2][l], cxx[3][l], cxx[4][l], cxx[5][l], lambda_th);
            }
    return s;
}


// ************************************************************************************
// given dark matter fields dm and Gaussian smoothing radius Rs, return the web classify result at point p0
// ************************************************************************************
//std::vector<int> environment(std::vector<Particle>& dm, double Rs, std::vector<Particle>& p0)
//{
//    force_base_type(0,1);
//    force_kernel_type(2);
//    auto begin = std::chrono::steady_clock::now();
//
//    auto sc = sfc_r2c(sfc(dm),true);
//    auto w = wft(2, 0);
//    auto cxx = tidal_tensor(sc, w);
//    auto env = web_classify(cxx,p0);
//    
//    fftw_free(sc);
//    delete[] w;
//
//    auto end = std::chrono::steady_clock::now();
//    std::cout << "Time difference EnvironmentClassify  = "
//    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
//    << "[ms]" << std::endl;
//
//    return env;
//}


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



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// we first split halo cataloge to four environmental sub_catalogues then in each sub-cata 
// further split according halo mass in envi major sequences, halo mass split is an independent parameter
std::vector<std::vector<Particle>*> halo_envi_mass_multi_split(std::vector<int>& envi, std::vector<Particle>& hl, int nbin)
{
    // ------envi split first--------
    auto evpts = halo_envi_match_and_split(envi, hl);

    // -----nodes of mass split-----
    std::vector<double> vecmass(hl.size());
    #pragma omp parallel for
    for(size_t i = 0; i < hl.size(); ++i) vecmass[i] = hl[i].weight;
    auto node = nodes_of_proto_sort(vecmass, nbin);

    // --------sub split-------
    std::vector<std::vector<Particle>*> cata;
    for(auto x : evpts){
        auto tmp = mass_classify_and_push_back(node,*x,nbin);
        for(auto y : tmp) cata.push_back(y);
    }

    return cata;
}

// unlike multi_split, we just concatenate splited sub-cataloges
std::vector<std::vector<Particle>*> halo_envi_mass_concatenate_split(std::vector<int>& envi, std::vector<Particle>& hl, int nbin)
{
    auto vpts = halo_envi_match_and_split(envi, hl);
    auto hvpts = halo_mass_split(hl,nbin);

    for(auto x : hvpts) vpts.push_back(x);

    return vpts;
}

// ************************************************************************************
// nbin is the number of catalogue splited, return as vector of catalog
// ************************************************************************************
std::vector<std::vector<Particle>*> halo_mass_split(std::vector<Particle>& hl, int nbin)
{
    std::vector<double> vecmass(hl.size());
    #pragma omp parallel for
    for(size_t i = 0; i < hl.size(); ++i) vecmass[i] = hl[i].weight;
    auto node = nodes_of_proto_sort(vecmass, nbin);

    return mass_classify_and_push_back(node,hl,nbin);
}

// given nodes of splitting halo cataloge, classify each particle with nodes and then
// push back to desire sub-catalog. parameter "node" returns from nodes_of_proto_sort()
std::vector<std::vector<Particle>*> mass_classify_and_push_back(std::vector<double>& node, std::vector<Particle>& hl, int nbin)
{
    std::vector<std::vector<Particle>*> cata;
    for(int i = 0; i < nbin; ++i) 
        cata.push_back(new std::vector<Particle>);

    for(auto x : hl){
        cata[classify_index(node,x.weight)]->push_back(x);
    }

    for(auto x : cata) if(x->size() - hl.size()/nbin > nbin) {
        std::cout << "[func: SPLIT] !Warning, some nodes include multiple identical items\n";
        break;
    }
    for(auto x : cata) print_min_max_and_size(*x);

    return cata;
}

void print_min_max_and_size(std::vector<Particle>& hl){
    if(hl.size() != 0){
        double sum{0};
        double min{hl[0].weight},max{hl[0].weight};
        for(auto x : hl){
            sum += x.weight;
            if(x.weight > max) max = x.weight;
            else if(x.weight < min) min = x.weight;
        }
        std::cout << "size: " << hl.size() << ", min: " << min << ", max: " << max << ", ave: " << sum/hl.size() << std::endl; 
    }
    else 
        std::cout << "!Empty vector\n";
}
void print_min_max_and_size_double(std::vector<double>& vec){
    if(vec.size() != 0){
        double sum{0};
        double min{vec[0]},max{vec[0]};
        for(auto x : vec){
            sum += x;
            if(x > max) max = x;
            else if(x < min) min = x;
        }
        std::cout << "size: " << vec.size() << ", min: " << min << ", max: " << max << ", ave: " << sum/vec.size() << std::endl; 
    }
    else 
        std::cout << "!Empty vector\n";
}
// which vector should trial been push back
//int classify_index(std::vector<double>& node, double trial){
//    int index{0};
//    for(int i = 0; i < node.size(); ++i){
//        if(trial > node[i]) ++index;
//        else break;
//    }
//    return index;
//}




// ************************************************************************************
// return the nbin fraction node points of a double vector in ascending order
// ************************************************************************************
std::vector<double> nodes_of_proto_sort(std::vector<double>& vec, int nbin)
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
std::vector<size_t> limited_sort(std::vector<double>& vec)
{
    std::vector<double> tmp;
    for(auto x : vec) tmp.push_back(x);
    std::cout << "[func: limited_sort] node size: " << tmp.size() << std::endl;
    std::vector<size_t> sortedID;
    size_t max_id = maximum_index(tmp);
    const double MAX{tmp[max_id]};
    for(size_t i = 0; i < tmp.size() - 1; ++i){
        size_t id = minimum_index(tmp);
        tmp[id] = MAX;
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
// return the optimal cross-correlation coefficients r and the corresponding weight vector 
// as {r_max,eigen_vector},(By solving A*x=lambda*C*x) which has size of n+1, input parameter 
// "cov" is returned by "covar_of_data_vector()", which stores the covariance of dm and n splited halo catalogues.  
// ************************************************************************************
std::vector<double> optimal_weight_solver(std::vector<double> cov, int n, bool PRINT)
{
    const int symsize = n*(n+1)/2;

    double A[symsize];
    double C[symsize];
    double Z[n*n];
    double lambda[n];
    //auto result = new double[n+1]{0};
    std::vector<double> result(n+1);

    //----initialize A and C-----
    int l = 0;
    for(int i = 0; i < n; ++i)
        for(int j = i; j < n; ++j)
        {
            A[l] = cov[i+1] * cov[j+1] / cov[0];
            ++l;
        }
    for(int i = 0; i < symsize; ++i) C[i] = cov[i+n+1];

    //----solving eigen value lambda and eigen vector Z-----
    auto info = LAPACKE_dspgv(LAPACK_ROW_MAJOR,1,'V','U',n,A,C,lambda,Z,n);

    //---print---
    if(info == 0 && PRINT){
        std::cout << "Eigen value: \n";
        for(int i = 0; i < n; ++i) std::cout << lambda[i] << " ";std::cout <<"\n";
        std::cout << "Eigen vector (in column): \n";
        double sum[n]{0};
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j) sum[j] += Z[i*n+j];
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j)
                std::cout << Z[i*n+j]/sum[j] << ", ";
            std::cout << "\n";
        }
    }

    //---if everything goes well, normalize eigenvector and return---
    if(info == 0 && lambda[n-1] <= 1.){
        result[0] = sqrt(lambda[n-1]);
        double sum{0};
        for(int i = 0; i < n; ++i) sum += Z[i*n+n-1];
        for(int i = 0; i < n; ++i) result[i+1] = Z[i*n+n-1] / sum;
    }

    std::cout << "----------------------------------------------------------\n";
    std::cout << "[Func: optimal_weight_solver] solution: \n" << "MAX(r) = " << result[0] << "\n";
    std::cout << "WEIGHT = ("; for(int i = 1; i < result.size()-1; ++i) std::cout << result[i] << ", ";
    std::cout << result[result.size()-1] << ")^T\n";
    std::cout << "----------------------------------------------------------\n";

    return result;
}

// variable {hlsize} for envi size matching
std::vector<int> envi_vector_readin(std::string ifn, size_t hlsize)
{
    std::ifstream ifs {ifn};
    if(!ifs){std::cout << "reading " + ifn + " with error, Abort"; std::terminate();}

    std::vector<int> envi;
    int temp{0};
    char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hlsize == envi.size()) 
        std::cout << "halo size matched! continue\n"; 
    else std::terminate();

    return envi;
}

std::vector<int> envi_with_Mcut(std::string ifn, double Mcut, std::vector<Particle>& hl)
{
    std::ifstream ifs {ifn};
    if(!ifs){std::cout << "reading " + ifn + " with error, Abort"; std::terminate();}

    std::vector<int> envi, enviMcut;
    int temp{0};
    char comma{0};
    while(ifs >> temp >> comma) envi.push_back(temp);
    if(hl.size() == envi.size()) 
        std::cout << "halo size matched! continue\n"; 
    else std::terminate();

    for(int i = 0; i < envi.size(); ++i){
        if(hl[i].weight > Mcut) enviMcut.push_back(envi[i]);
    }
    std::vector<int>().swap(envi);

    return enviMcut;
}

// ************************************************************************************
// this function return the environmental split halo catalogue and dm as {dm,vd,st,fl,kt}
// specialized for optimal weight solver 
// ************************************************************************************
std::vector<std::vector<Particle>*> halo_envi_match_and_split(std::vector<int>& envi, std::vector<Particle>& hl)
{
    auto vd = new std::vector<Particle>;
    auto st = new std::vector<Particle>;
    auto fl = new std::vector<Particle>;
    auto kt = new std::vector<Particle>;
    for(size_t i = 0; i < envi.size(); ++i){
        if(envi[i] == 0) vd->push_back(hl[i]);
        else if(envi[i] == 1) st->push_back(hl[i]);
        else if(envi[i] == 2) fl->push_back(hl[i]);
        else if(envi[i] == 3) kt->push_back(hl[i]);
    }

    //std::vector<int>().swap(envi);

    if(hl.size() != (vd->size() + st->size() + fl->size() + kt->size())) 
        std::cout << "Warning! halo environment subset size not matched\n";
    std::vector<std::vector<Particle>*> result{vd,st,fl,kt};

    for(auto x : result) print_min_max_and_size(*x);
    
    return result;
}

// ************************************************************************************************************
// return the fourier of scaling function coefficients of optimal reconstructed halo catalogues, vpts is the
// vector of dm and splited halo cataloges, in envi-split case: {dm,vd,st,fl,kt}. R specify the smmothing scale
// which decide the reconstruct coeefficient of each halo component (weight vector), after solving weight
// vector we then reconstruct the halo fields with optimal weight and return as fourier of sfc coefficients
// ************************************************************************************************************
fftw_complex* optimal_reconstruct(fftw_complex* sc_dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT)
{
    // ------covariance of each component------
    std::vector<double> cov;

    auto wpk = window_Pk(R,0);

    std::vector<fftw_complex*> vec_sc; vec_sc.push_back(sc_dm);

    // -------size check before eigen solver-----
    bool EmptySize {false};
    const int ThdSize {100}; 

    for(auto x : vpts) 
        if (x->size() < ThdSize) 
            EmptySize = true;

    if(EmptySize){
        std::cout << "!EmptySize, some elements will be removed\n";
        std::cout << "---+bf\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
        for(int i = 0; i < vpts.size(); ++i) {
            if(vpts[i]->size() < ThdSize){
                delete vpts[i];
                vpts.erase(vpts.begin()+i);
                --i;
            }
        }
        std::cout << "---+af\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
    }

    // -------covariance array-------
    for(auto x : vpts) vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
        {   
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk,true));
        }

    // ------solving optimal weight vector------
    size_t dim{vpts.size()}; 
    auto solve =  optimal_weight_solver(cov,dim,PRINT);

    // -----------------reconstruct--------------
    std::vector<double> weight, norm;
    for(int i = 1; i < solve.size(); ++i) weight.push_back(solve[i]);
    for(auto x : vpts){
        double sum = 0;
        //#pragma omp parallel reduction (+:sum)
        for(auto pt : *x) sum += pt.weight;
        norm.push_back(sum);
    }

    double total{0};
    for(auto x : norm) total += x;
    for(int i = 0; i < weight.size(); ++i) weight[i] /= norm[i] / total;

    #pragma omp parallel for
    for(size_t l = 0; l < GridLen * GridLen * (GridLen/2 + 1); ++l) {
            vec_sc[1][l][0] *= weight[0];
            vec_sc[1][l][1] *= weight[0];
        }
    
    for(int i = 2; i < vec_sc.size(); ++i){
        #pragma omp parallel for
        for(size_t l = 0; l < GridLen * GridLen * (GridLen/2 + 1); ++l) {
            vec_sc[1][l][0] += vec_sc[i][l][0] * weight[i-1];
            vec_sc[1][l][1] += vec_sc[i][l][1] * weight[i-1];
        }
        fftw_free(vec_sc[i]);
    }

    return vec_sc[1];
}


// vector of {sqrt(lambda),eigen_vector}
std::vector<double> optimal_solution(std::vector<Particle>& dm, std::vector<std::vector<Particle>*> vpts, double R, bool PRINT)
{
    // ------covariance of each component------
    std::vector<double> cov;

    auto wpk = window_Pk(R,0);

    std::vector<fftw_complex*> vec_sc; vec_sc.push_back(sfc_r2c(sfc(dm),true));

    // -------size check before eigen solver-----
    bool EmptySize {false};

    for(auto x : vpts) 
        if (x->size() < 100) 
            EmptySize = true;

    if(EmptySize){
        std::cout << "!EmptySize, some elements will be removed\n";
        std::cout << "---+bf\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
        for(int i = 0; i < vpts.size(); ++i) {
            if(vpts[i]->size() < 100){
                delete vpts[i];
                vpts.erase(vpts.begin()+i);
                --i;
            }
        }
        std::cout << "---+af\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
    }

    // -------covariance array-------
    for(auto x : vpts) vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
        {   
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk,true));
        }
    for(auto x : vec_sc) fftw_free(x);

    // ------solving optimal weight vector------
    size_t dim{vpts.size()}; 
    auto solve =  optimal_weight_solver(cov,dim,PRINT);

    return solve;
}

// wpk is return by window_pk()
std::vector<double> optimal_solution_lean(fftw_complex* sc_dm, std::vector<std::vector<Particle>*> vpts, double* wpk, bool PRINT)
{
    // ------covariance of each component------
    std::vector<double> cov;

    std::vector<fftw_complex*> vec_sc; vec_sc.push_back(sc_dm);

    // -------size check before eigen solver-----
    bool EmptySize {false};

    for(auto x : vpts) 
        if (x->size() < 100) 
            EmptySize = true;

    if(EmptySize){
        std::cout << "!EmptySize, some elements will be removed\n";
        std::cout << "---+bf\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
        for(int i = 0; i < vpts.size(); ++i) {
            if(vpts[i]->size() < 100){
                delete vpts[i];
                vpts.erase(vpts.begin()+i);
                --i;
            }
        }
        std::cout << "---+af\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
    }

    // -------covariance array-------
    for(auto x : vpts) vec_sc.push_back(sfc_r2c(sfc(*x),true));

    for(int i = 0; i < vec_sc.size(); ++i)
        for(int j = i; j < vec_sc.size(); ++j)
        {   
            cov.push_back(covar_CombinewithKernel(densityCovarianceArray(vec_sc[i],vec_sc[j]),wpk,true));
        }
    for(int i = 1; i < vec_sc.size(); ++i) fftw_free(vec_sc[i]);

    // ------solving optimal weight vector------
    size_t dim{vpts.size()}; 
    auto solve =  optimal_weight_solver(cov,dim,PRINT);

    return solve;
}
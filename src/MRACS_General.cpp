// =======================================================
// This cpp source file is part of MRACS project
// 
// =======================================================
#include"MRACS_General.h"


std::vector<double> log_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT)
{
    if(Rmin <= 0) {std::cout << "input error, Rmin should be positive number\n";std::terminate();}
    std::vector<double> r(Npt);
    for(int i = 0; i < Npt; ++i) r[i] = Rmin * pow(Rmax/Rmin, static_cast<double>(i)/Npt);
    if(ENDPOINT) r.push_back(Rmax);

    return r;
}

std::vector<double> linear_scale_generator(double Rmin, double Rmax, int Npt, bool ENDPOINT)
{
    std::vector<double> r(Npt);
    for(int i = 0; i < Npt; ++i) r[i] = Rmin + i * (Rmax - Rmin) / Npt;
    if(ENDPOINT) r.push_back(Rmax);

    return r;
}

// generate n random particle locate in box (0,boxsize)^3 
std::vector<Particle> default_random_particle(double boxsize, size_t n)
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> v(0,boxsize);

    std::vector<Particle> p0;
    for(size_t i = 0; i < n; ++i) 
        p0.push_back({v(e), v(e), v(e), 1.});
    
    return p0;
}

// generate n random particle locate in box (Lx,Ly,Lz) 
std::vector<Particle> default_random_particle(double Lx, double Ly, double Lz, size_t n)
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> vx(0,Lx);
    std::uniform_real_distribution<double> vy(0,Ly);
    std::uniform_real_distribution<double> vz(0,Lz);

    std::vector<Particle> p0;
    for(size_t i = 0; i < n; ++i) 
        p0.push_back({vx(e), vy(e), vz(e), 1.});
    
    return p0;
}

// x is the number of particles per side, i.e there are x^3 particles in all
// L is the box size and w is the safe band width, i.e random point locate in (w,L-w)^3
std::vector<Particle> generate_random_particle(int x, double L, double w)
{
    double diff;
    clock_t begin, end;

    std::vector<Particle> p;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0,1);
    const int64_t nps = x;
    const int64_t NPt = nps * nps * nps;
    const double boxL = L;
    const double safeband = w;

    begin = clock();
    for(int i = 0; i < nps; ++i)
        for(int j = 0; j < nps; ++j)
            for(int k = 0; k < nps; ++k)
                p.push_back({safeband + (i + u(e)) * (boxL - 2*safeband) / nps,
                              safeband + (j + u(e)) * (boxL - 2*safeband) / nps,
                              safeband + (k + u(e)) * (boxL - 2*safeband) / nps, 1.});
    end = clock();
    diff = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "generating time for " << p.size() << " random points:  " << diff << "s" << std::endl;
    
    return p;
}


//calculate index of box that must be inside the sphere R, and that might be intersect with sphere boundray
void fill_index_set(const double R, std::vector<Index>& inner_index, std::vector<Index>& cross_index)
{
    const int M = R + 1;
    for(int i = -M; i <= M; ++i)
        for(int j = -M; j <= M; ++j)
            for(int k = -M; k <= M; ++k)
            {
                if(sqrt(i * i + j * j + k * k) <= R - sqrt(3.))
                {
                    inner_index.push_back(Index{i, j, k});
                }
                else if(sqrt(i * i + j * j + k * k) < R + sqrt(3.))
                {
                    cross_index.push_back(Index{i, j, k});
                }
            }
}


//particle periodic in box size L, calculate the number of particles in sphere center at p0 with radius R
double* count_in_sphere(const double R, std::vector<Particle>& p, std::vector<Particle>& p0)
{
    std::vector<Index> inner_index;
    std::vector<Index> cross_index;
    std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();

    fill_index_set(R/SimBoxL, inner_index, cross_index);

    auto count = new double[p0.size()];
    double temp;

    #pragma omp parallel for reduction(+:temp)

    for(size_t n = 0; n < p0.size(); ++n)
    {
        temp = 0;
        for(size_t i = 0; i < p.size(); ++i)
            for(size_t m = 0; m < cross_index.size(); ++m)
            {
                double xx = p[i].x + cross_index[m].i * SimBoxL - p0[n].x;
                double yy = p[i].y + cross_index[m].j * SimBoxL - p0[n].y;
                double zz = p[i].z + cross_index[m].k * SimBoxL - p0[n].z;
                if((abs(xx) < R) && (abs(yy) < R) && (abs(zz) < R))
                    if(xx*xx + yy*yy + zz*zz < R*R)
                        ++temp;
            }
        count[n]= temp + inner_index.size() * p.size();
    }

    std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 count    = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return count;
}


// safe band is needed, i.e, no periodic boundary condictions added in
double* count_in_cylinder(double R, double H, std::vector<Particle>& p, std::vector<Particle>& p0)
{
    std::chrono::steady_clock::time_point begin4 = std::chrono::steady_clock::now();
    auto count = new double[p0.size()];
    double temp;

    #pragma omp parallel for reduction(+:temp)
    for(size_t n = 0; n < p0.size(); ++n)
    {
        temp = 0;
        for(size_t i = 0; i < p.size(); ++i)
        {
            double xx = p[i].x - p0[n].x;
            double yy = p[i].y - p0[n].y;
            double zz = p[i].z - p0[n].z;
            if((abs(xx) < R) && (abs(yy) < R) && (abs(zz) < H/2))
                if(xx*xx + yy*yy  < R*R)
                    ++temp;
        }
            
        count[n]= temp;
    }

    std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
    std::cout << "Time difference 4 count_cyl = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - begin4).count()
    << "[ms]" << std::endl;

    return count;
}

// sampling the window fourier from [1e-4,10] with delta_k = 1e-4
// combine with Pk_theory value to calculate the theoretical second order statistics
// kernel can be "shell" "Gaussian","sphere","GDW","GLS"
std::vector<double> win_theory(std::string kernel, double R, double theta)
{
    const double k_min = 1e-4; //radian frequency
    const double k_max = 10;
    const double dk = 1e-4;
    const int klen = 1e5;

    std::vector<double> win(klen);
    if(kernel == "shell") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = sin(phase)/(phase);
        }
    }
    else if(kernel == "sphere") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = 3*(sin(phase)-phase*cos(phase))/(pow(phase,3));
        }
    }
    else if(kernel == "Gaussian") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = pow(1/M_E,phase*phase/2);
        }
    }
    else if(kernel == "GDW") {
        const double norm = pow(2,7./4)/sqrt(15)*pow(TWOPI,3./4)*pow(R,3./2);
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase = (i+1)*dk*R;
            win[i] = norm * pow(phase,2) * pow(1/M_E,phase*phase/2);
        }
    }
    else if(kernel == "ThickShell") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase1 = (i+1)*dk*R;
            double phase2 = (i+1)*dk*theta;
            win[i] = 3*(sin(phase2)-phase2*cos(phase2)-sin(phase1)+phase1*cos(phase1))/(pow(phase2,3)-pow(phase1,3));
        }
    }
    else if(kernel == "GLS") {
        #pragma omp parallel for
        for(int i = 0; i < klen; ++i){
            double phase1 = (i+1)*dk*R;
            double phase2 = (i+1)*dk*theta;
            win[i] = (theta*theta * cos(phase1) + R*R * sin(phase1)/phase1)
                     /(R*R + theta*theta) * pow(1/M_E,phase2 * phase2 / 2);
        }
    }
    return win;
}